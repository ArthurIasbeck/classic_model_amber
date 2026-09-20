from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt

import ross as rs
import control as ct
from scipy.linalg import eigh, block_diag
from scipy.optimize import least_squares

from ross import MagneticBearingElement


def get_rotor():
    ## Shaft material creation ##
    steel = rs.Material(name="steel", rho=7850, E=1.9e11, Poisson=0.30)
    steel_m12 = rs.Material(name="steel", rho=7700, E=2e11, Poisson=0.31, color="red")

    ## Shaft elements ##
    Li = [
        0.0,
        0.012,
        0.032,
        0.052,
        0.072,
        0.092,
        0.112,
        0.1208,
        0.12724,
        0.13475,
        0.14049,
        0.14689,
        0.15299,
        0.159170,
        0.16535,
        0.180350,
        0.1905,
        0.2063,
        0.2221,
        0.2379,
        0.2537,
        0.2695,
        0.2853,
        0.3011,
        0.3169,
        0.3243,
        0.3363,
        0.358,
        0.364,
        0.3705,
        0.3825,
        0.3986,
        0.4147,
        0.4308,
        0.4469,
        0.4630,
        0.4791,
        0.4952,
        0.5113,
        0.5274,
        0.5356,
        0.5457,
        0.5607,
        0.5669,
        0.5731,
        0.5792,
        0.5856,
        0.5913,
        0.5989,
        0.6053,
        0.6141,
        0.6341,
        0.6461,
    ]
    # Shaft discretization - node positions ## CG: 0.3243
    Li = [round(i, 4) for i in Li]  # Rounding decimal places
    L = [Li[i + 1] - Li[i] for i in range(len(Li) - 1)]  # Element size (e = n-1)
    i_d = [0.0 for _ in L]  # Internal diameter
    o_d1 = [0.0 for _ in L]  # External diameter

    # Adjustments for external diameter
    o_d1[0] = 6.35
    o_d1[1:5] = [32 for _ in range(4)]
    o_d1[5:14] = [34.8 for _ in range(9)]
    o_d1[14:16] = [49.9 for _ in range(2)]
    o_d1[16:27] = [19.05 for _ in range(11)]
    o_d1[27:29] = [54 for _ in range(2)]
    o_d1[29:40] = [19.05 for _ in range(12)]
    o_d1[40:42] = [49.9 for _ in range(2)]
    o_d1[42:51] = [34.8 for _ in range(9)]
    o_d1[51] = 6.35
    o_d = [i * 1e-3 for i in o_d1]  # Conversion to meters

    shaft_elements = [
        rs.ShaftElement(
            L=l,
            idl=idl,
            odl=odl,
            material=steel,
            shear_effects=True,
            rotary_inertia=True,
            gyroscopic=True,
            alpha=4,
            beta=1e-3,
        )
        for l, idl, odl in zip(L, i_d, o_d)
    ]

    # Disk elements ##
    # n_list = [27, 28, 29]  # Central disk positioning
    n_list = [28]  # Central disk positioning
    n_list_2 = [
        6,
        7,
        8,
        9,
        10,
        11,
        12,
        13,
        43,
        44,
        45,
        46,
        47,
        48,
        49,
        50,
    ]  # Disks representing lamination
    width = [0.004, 0.007, 0.014]  # Central disk width

    width_2 = [
        0.0088,
        0.0064,
        0.0075,
        0.0057,
        0.0064,
        0.0061,
        0.0062,
        0.0062,
        0.0062,
        0.0062,
        0.0061,
        0.0064,
        0.0057,
        0.0075,
        0.0064,
        0.0088,
    ]  # Lamination disk width

    i_disc_1 = [0.054, 0.054, 0.054]  # Central disk internal diameter
    i_disc_2 = [0.0348] * 16  # Lamination disk internal diameter
    o_disc = [0.1200] * 3  # Central disk external diameter
    o_disc_2 = [0.0498] * 16  # Lamination disk external diameter

    disk_elements_1 = [
        rs.DiskElement.from_geometry(n=n, material=steel, width=m, i_d=Id, o_d=Od)
        for n, m, Id, Od in zip(n_list, width, i_disc_1, o_disc)
    ]

    disk_elements_2 = [
        rs.DiskElement.from_geometry(n=n, material=steel_m12, width=m, i_d=Id, o_d=Od)
        for n, m, Id, Od in zip(n_list_2, width_2, i_disc_2, o_disc_2)
    ]

    disk_elements = [*disk_elements_1, *disk_elements_2]

    # Bearing elements:
    s = ct.TransferFunction.s
    n_list = [12, 43]
    n = 138
    A = 470.3e-6
    i0 = 1.0
    s0 = 0.432e-3
    alpha = 0.392
    # c_13 = (
    #     1e6
    #     * 0.0062
    #     * (s + 46)
    #     / s
    #     * 0.0062
    #     * (400 / 77)
    #     * ((s + 77 * 2 * np.pi) / (s + 400 * 2 * np.pi))
    #     * (409 / 124)
    #     * ((s + 124 * 2 * np.pi) / (s + 409 * 2 * np.pi))
    #     * ct.tf(
    #         [1, 371.964570185032, 5404595.37003653],
    #         [1, 1301.87599564761, 3095107.94018162],
    #     )
    #     * ct.tf(
    #         [1, 282.743338823081, 22206609.9024511],
    #         [1, 1337.06183336782, 17458343.225087],
    #     )
    #     * 0.4
    # )

    c_24 = (
        1e6
        * 0.0046
        * 0.0046
        * (s + 35)
        / s
        * (75 / 25)
        * ((s + 25 * 2 * np.pi) / (s + 75 * 2 * np.pi))
        * (1690 / 260)
        * ((s + 260 * 2 * np.pi) / (s + 1690 * 2 * np.pi))
        * ct.tf(
            [1, 1468.38040628787, 3206634.46991393],
            [1, 931.168062524015, 3206634.46991393],
        )
        * ct.tf(
            [1, 1357.16802635079, 6316546.81669719],
            [1, 1240.30077963725, 8720782.44880256],
        )
        * 2
    )
    k_amp = 1.0
    k_sense = 1.0
    bearing_elements = [
        # O mancal será comentado no modelo porque será representado por meio da aplicação de forças eletromagnéticas.
        # MagneticBearingElement(
        #     n=n_list[0],
        #     g0=s0,
        #     i0=i0,
        #     ag=A,
        #     nw=n,
        #     alpha=alpha,
        #     controller_transfer_function=c_13,
        #     k_amp=k_amp,
        #     k_sense=k_sense,
        # ),
        # rs.MagneticBearingElement(
        #     n=n_list[1],
        #     g0=s0,
        #     i0=i0,
        #     ag=A,
        #     nw=n,
        #     alpha=alpha,
        #     controller_transfer_function=c_24,
        #     k_amp=k_amp,
        #     k_sense=k_sense,
        # ),
    ]

    mma = rs.MagneticBearingElement(
            n=n_list[1],
            g0=s0,
            i0=i0,
            ag=A,
            nw=n,
            alpha=alpha,
            controller_transfer_function=c_24,
            k_amp=k_amp,
            k_sense=k_sense,
        )
    print(f"ki = {mma.ki}")
    print(f"ks = {mma.ks}")

    ## Rotor assembly - 6dof ##
    # Rotor construction:
    rotor = rs.Rotor(
        shaft_elements=shaft_elements,
        disk_elements=disk_elements,
        bearing_elements=bearing_elements,
    )

    return rotor


class BuildMimoModelFromFrequencyResponse:
    def __init__(self):
        self.G = None
        self.G_w13_response = None
        self.G_v13_response = None
        self.angular_frequencies = None
        self.rotor = get_rotor()
        self.x = None
        self.n_u = None
        self.n_x = None
        self.node = None
        self.n_amb = None
        self.n_dof = None
        self.k_x = None
        self.k_s = None
        self.k_i = None
        self.original_model_v = None
        self.original_model_w = None
        self.G_v = None
        self.G_w = None
        self.speed = 0
        self.amb = None
        self.num_modes = 20

        self.Phi = None
        self.M_m = None
        self.C_m = None
        self.K_m = None

        self.data_dir = "../data"

    def process_rotor(self):
        self.n_dof = self.rotor.ndof
        self.n_dof = self.rotor.ndof
        self.n_x = 2 * self.n_dof
        self.n_u = self.n_dof
        self.x = np.zeros((self.n_x, 1))  # TODO: Verificar se é necessário

    def setup_modal_domain(self):
        self.process_rotor()

        modal_reduction = True

        M = self.rotor.M(self.speed)
        C = self.rotor.C(self.speed)
        K = self.rotor.K(self.speed)

        eigenvalues, eigenvectors = eigh(K, M)

        if self.num_modes == -1:
            modal_reduction = False
            self.num_modes = len(eigenvalues)

        selected_eigenvalues = eigenvalues[: self.num_modes]

        omega_rad_s = np.sqrt(np.abs(selected_eigenvalues))
        freq_hz = omega_rad_s / (2 * np.pi)
        speed_rpm = freq_hz * 60

        if modal_reduction:
            print("\n" + 120 * "=")
            print(f"- Frequencies of the first {self.num_modes} modes")
            for i in range(self.num_modes):
                print(f"Mode {i + 1}: {freq_hz[i]:.6f} Hz | {speed_rpm[i]:.2f} RPM")
            print(120 * "=")

        self.Phi = eigenvectors[:, : self.num_modes]
        self.M_m = self.Phi.T @ M @ self.Phi
        self.C_m = self.Phi.T @ C @ self.Phi
        self.K_m = self.Phi.T @ K @ self.Phi

    def _build_model_response(self, k_x, k_i, show_pole_zero_map=False):
        """Build the reduced model and return its complex frequency response."""
        theta = np.pi / 4  # Bearing orientation angle (45 degrees)
        n_controllers = 2

        Phi = self.Phi

        phi = np.zeros((self.n_dof, n_controllers))
        phi[12 * 6 + 0, 0] = 1  # x
        phi[12 * 6 + 1, 1] = 1  # y

        c_theta = np.cos(theta)
        s_theta = np.sin(theta)
        R = np.array([[c_theta, s_theta], [-s_theta, c_theta]])

        K_x = self._gain_matrix(k_x, "k_x")
        K_i = self._gain_matrix(k_i, "k_i")

        N = self.n_dof
        n_c = n_controllers
        m = self.num_modes

        zeros = lambda rows, columns: np.zeros((rows, columns))

        H = np.block([np.eye(m), zeros(m, m)])
        phi_t = np.transpose(phi)

        M_m_I = np.linalg.inv(self.M_m)
        Phi_T = np.transpose(self.Phi)
        A = np.block([[zeros(m, m), np.eye(m)], [-M_m_I @ self.K_m, -M_m_I @ self.C_m]])
        B = np.block([[zeros(m, N)], [M_m_I @ Phi_T]])

        A_star = A + B @ phi @ K_x @ R @ phi_t @ Phi @ H
        B_star = B @ phi @ K_i
        C_star = 1e6 * R @ phi_t @ Phi @ H
        D_star = zeros(n_c, n_c)

        model = ct.ss(A_star, B_star, C_star, D_star)

        if show_pole_zero_map:
            model_for_pole_zero_map = ct.minreal(model, tol=1e-3, verbose=False)
            poles = ct.poles(model_for_pole_zero_map)
            zeros = ct.zeros(model_for_pole_zero_map)

            fig, ax = plt.subplots(figsize=(6, 5), dpi=180)
            ax.plot(poles.real, poles.imag, "x", ms=8, mew=1.5, label="Poles")
            ax.plot(
                zeros.real,
                zeros.imag,
                "o",
                ms=7,
                mfc="none",
                mew=1.5,
                label="Zeros",
            )
            ax.axhline(0, color="black", linewidth=0.8)
            ax.axvline(0, color="black", linewidth=0.8)
            ax.set_xlabel("Real axis")
            ax.set_ylabel("Imaginary axis")
            ax.set_title("Pole-zero map of the reduced model")
            ax.grid(True)
            ax.legend()

            fig.tight_layout()

        model_mag, model_phase, _ = ct.frequency_response(
            model, self.angular_frequencies
        )
        return model_mag * np.exp(1j * model_phase)

    @staticmethod
    def _gain_matrix(gains, name):
        gains = np.asarray(gains, dtype=float)
        if gains.shape == (2,):
            return np.diag(gains)
        if gains.shape == (2, 2):
            return gains
        raise ValueError(f"{name} must contain two gains or be a 2x2 matrix.")

    def _experimental_responses(self):
        if (
            self.angular_frequencies is None
            or self.G_v13_response is None
            or self.G_w13_response is None
        ):
            raise ValueError(
                "Experimental frequency responses must be loaded before model fitting."
            )

        omega = np.asarray(self.angular_frequencies).squeeze()
        responses = (
            np.asarray(self.G_v13_response).squeeze(),
            np.asarray(self.G_w13_response).squeeze(),
        )
        if omega.ndim != 1 or any(response.ndim != 1 for response in responses):
            raise ValueError("Frequencies and experimental responses must be 1-D arrays.")
        if any(response.size != omega.size for response in responses):
            raise ValueError(
                "Frequencies and experimental responses must have the same length."
            )
        return omega, responses

    def optimize_gains(
        self,
        initial_guess=(129012.78635407916, 129012.78635407916, 55.733523704962195, 55.733523704962195),
        # bounds=((-1e6, -1e6, 0.0, 0.0), (0.0, 0.0, 1e4, 1e4)),
        magnitude_weight=1.0,
        phase_weight=1.0,
        max_nfev=300,
        verbose=2,
    ):
        _, experimental = self._experimental_responses()
        self.setup_modal_domain()

        initial_guess = np.asarray(initial_guess, dtype=float)
        # lower, upper = (np.asarray(bound, dtype=float) for bound in bounds)
        # if initial_guess.shape != (4,):
        #     raise ValueError("initial_guess must contain four diagonal gains.")
        # if lower.shape != (4,) or upper.shape != (4,):
        #     raise ValueError("bounds must contain four lower and four upper limits.")
        # if np.any(initial_guess < lower) or np.any(initial_guess > upper):
        #     raise ValueError("initial_guess must be inside bounds.")
        # if magnitude_weight < 0 or phase_weight < 0:
        #     raise ValueError("Residual weights must be non-negative.")
        # if magnitude_weight == 0 and phase_weight == 0:
        #     raise ValueError("At least one residual weight must be positive.")

        experimental_magnitude = np.abs(np.asarray(experimental))
        experimental_phase = np.angle(np.asarray(experimental))
        epsilon = np.finfo(float).tiny

        def residual(gains):
            model_response = self._build_model_response(gains[:2], gains[2:])
            model_response = (model_response[0, 0, :], model_response[1, 1, :])
            model_magnitude = np.abs(np.asarray(model_response))
            model_phase = np.angle(np.asarray(model_response))

            magnitude_residual = 20 * np.log10(
                np.maximum(model_magnitude, epsilon)
            ) - 20 * np.log10(np.maximum(experimental_magnitude, epsilon))
            phase_residual = np.rad2deg(
                np.angle(np.exp(1j * (model_phase - experimental_phase)))
            )
            return np.concatenate(
                [
                    magnitude_weight * magnitude_residual.ravel(),
                    phase_weight * phase_residual.ravel(),
                ]
            )

        result = least_squares(
            residual,
            initial_guess,
            # bounds=(lower, upper),
            max_nfev=max_nfev,
            x_scale="jac",
            verbose=verbose,
        )
        self.k_x = np.diag(result.x[:2])
        self.k_s = self.k_x
        self.k_i = np.diag(result.x[2:])
        if verbose:
            print("\nOptimized K_x:")
            print(self.k_x)
            print("Optimized K_i:")
            print(self.k_i)
        return result, self.k_x, self.k_i

    def build_model(self, k_x=None, k_i=None):
        omega = self.angular_frequencies
        self._experimental_responses()
        self.setup_modal_domain()

        if k_x is None:
            k_x = (129012.78635407916, 129012.78635407916)
        if k_i is None:
            k_i = (55.733523704962195, 55.733523704962195)
        self.k_x = self._gain_matrix(k_x, "k_x")
        self.k_s = self.k_x
        self.k_i = self._gain_matrix(k_i, "k_i")
        model_response = self._build_model_response(
            k_x, k_i, show_pole_zero_map=True
        )

        model_v13 = model_response[0, 0, :]
        model_phase_v13 = np.angle(model_v13)

        model_w13 = model_response[1, 1, :]
        model_phase_w13 = np.angle(model_w13)

        exp_mag_v13 = np.abs(self.G_v13_response)
        exp_phase_v13 = np.angle(self.G_v13_response)

        exp_mag_w13 = np.abs(self.G_w13_response)
        exp_phase_w13 = np.angle(self.G_w13_response)

        fig = plt.figure(figsize=(6, 5), dpi=180)
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, 20 * np.log10(np.abs(model_v13)), label="Model")
        plt.semilogx(omega, 20 * np.log10(exp_mag_v13), label="Experimental")
        plt.legend()
        plt.grid()
        plt.ylabel("Magnitude [dB]")

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, np.rad2deg(np.unwrap(model_phase_v13)))
        plt.semilogx(omega, np.rad2deg(np.unwrap(exp_phase_v13)))
        plt.grid()
        plt.xlabel("Angular frequency [rad/s]")
        plt.ylabel("Phase [deg]")

        fig.suptitle(
            rf"Open Loop Frequency Response $\rightarrow$ V13"
        )

        fig = plt.figure(figsize=(6, 5), dpi=180)
        plt.subplot(2, 1, 1)
        plt.semilogx(omega, 20 * np.log10(np.abs(model_w13)), label="Model")
        plt.semilogx(omega, 20 * np.log10(exp_mag_w13), label="Experimental")
        plt.legend()
        plt.grid()
        plt.ylabel("Magnitude [dB]")

        plt.subplot(2, 1, 2)
        plt.semilogx(omega, np.rad2deg(np.unwrap(model_phase_w13)))
        plt.semilogx(omega, np.rad2deg(np.unwrap(exp_phase_w13)))
        plt.grid()
        plt.xlabel("Angular frequency [rad/s]")
        plt.ylabel("Phase [deg]")

        fig.suptitle(
            rf"Open Loop Frequency Response $\rightarrow$ W13"
        )

    def get_open_loop_frequency_responses(self):
        data_v13 = np.load(
            Path(self.data_dir) / f"T_v13_frequency_response_open_loop.npz"
        )
        data_w13 = np.load(
            Path(self.data_dir) / f"T_w13_frequency_response_open_loop.npz"
        )

        self.angular_frequencies = data_v13["angular_frequencies"]
        self.G_v13_response = data_v13["response"]
        self.G_w13_response = data_w13["response"]

        G = []
        for i in range(self.angular_frequencies.size):
            G.append(np.block([[self.G_v13_response[i], 0], [0, self.G_w13_response[i]]]))

        G = np.array(G)
        self.G = G

        output_path = Path(self.data_dir) / "open_loop_frequency_responses.txt"
        np.savetxt(
            output_path,
            np.column_stack(
                (
                    np.asarray(self.angular_frequencies),
                    np.real(self.G_v13_response),
                    np.imag(self.G_v13_response),
                    np.real(self.G_w13_response),
                    np.imag(self.G_w13_response),
                )
            ),
            fmt="%.18e",
        )


def main():
    # rotor = get_rotor()
    # figure = rotor.plot_rotor(nodes=999)
    # output_path = Path(__file__).resolve().parent.parent / "plots" / "rotor.png"
    # figure.write_image(output_path, scale=10)
    # figure.show()

    build_model = BuildMimoModelFromFrequencyResponse()
    build_model.get_open_loop_frequency_responses()
    build_model.build_model()

    # result, K_x, K_i = build_model.optimize_gains()
    # build_model.build_model(K_x, K_i)


if __name__ == "__main__":
    main()
    plt.show()
