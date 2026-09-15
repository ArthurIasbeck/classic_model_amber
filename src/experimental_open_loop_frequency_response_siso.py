from pathlib import Path

import control as ct
import numpy as np
from matplotlib import pyplot as plt


def get_controller():
    csiN1_W13_DE = 0.08
    fN1_W13_DE = 2.324778563656447e3

    csiD1_W13_DE = 0.37
    fD1_W13_DE = 1.759291886010284e3

    csiN2_W13_DE = 0.03
    fN2_W13_DE = 4.712388980384690e3

    csiD2_W13_DE = 0.16
    fD2_W13_DE = 4.178318229274425e3

    s = ct.tf("s")

    C = (
        0.0062
        * (1 + 46 / s)
        * (400 / 77)
        * ((s + 77 * 2 * np.pi) / (s + 400 * 2 * np.pi))
        * (409 / 124)
        * ((s + 124 * 2 * np.pi) / (s + 409 * 2 * np.pi))
        * ct.tf(
            [1, 2 * csiN1_W13_DE * fN1_W13_DE, fN1_W13_DE**2],
            [1, 2 * csiD1_W13_DE * fD1_W13_DE, fD1_W13_DE**2],
        )
        * ct.tf(
            [1, 2 * csiN2_W13_DE * fN2_W13_DE, fN2_W13_DE**2],
            [1, 2 * csiD2_W13_DE * fD2_W13_DE, fD2_W13_DE**2],
        )
        * 0.4
    )

    print(C)
    return C


class ExperimentalOpenLoopFrequencyResponseSiso:
    def __init__(self, data_name):
        self.H = None  # Open-loop frequency response
        self.data_name = data_name
        self.plots_dir = "../plots"
        self.data_dir = "../data"

        data = np.load(Path(self.data_dir) / f"{self.data_name}_frequency_response.npz")

        self.close_loop_fr = data["response"]
        self.angular_frequencies = data["angular_frequencies"]

        angular_freq_mask = (self.angular_frequencies > 4) & (
            self.angular_frequencies < 600
        )
        self.angular_frequencies = self.angular_frequencies[angular_freq_mask]
        self.close_loop_fr = self.close_loop_fr[angular_freq_mask]

    def compute_open_loop_freq_resp(self, plot=True):
        C = get_controller()
        omega = self.angular_frequencies
        controller_magnitude, controller_phase, _ = ct.frequency_response(C, omega)

        C_jw = controller_magnitude * np.exp(1j * controller_phase)
        T_jw = self.close_loop_fr
        G_jw = -T_jw / (C_jw * (1 + T_jw))

        magnitude = 20 * np.log10(np.abs(G_jw))
        phase = np.unwrap(np.angle(G_jw))
        self.H = G_jw

        np.savez(
            Path(self.data_dir) / f"{self.data_name}_frequency_response_open_loop.npz",
            angular_frequencies=self.angular_frequencies,
            response=self.H,
        )

        if plot:
            fig = plt.figure(figsize=(6, 5), dpi=180)
            plt.subplot(2, 1, 1)
            plt.semilogx(omega[omega < 600], 20 * np.log10(controller_magnitude))
            plt.grid()
            plt.ylabel("Magnitude [dB]")

            plt.subplot(2, 1, 2)
            plt.semilogx(omega, np.rad2deg(controller_phase))
            plt.grid()
            plt.xlabel("Angular frequency [rad/s]")
            plt.ylabel("Phase [deg]")

            fig.suptitle(f"Controller Frequency Response")
            plt.savefig(
                Path(self.plots_dir)
                / f"controller_frequency_response_{self.data_name}.svg",
                format="svg",
            )

            fig = plt.figure(figsize=(6, 5), dpi=180)
            plt.subplot(2, 1, 1)
            plt.semilogx(omega, magnitude)
            plt.grid()
            plt.ylabel("Magnitude [dB]")

            plt.subplot(2, 1, 2)
            plt.semilogx(omega, np.rad2deg(phase))
            plt.grid()
            plt.xlabel("Angular frequency [rad/s]")
            plt.ylabel("Phase [deg]")

            fig.suptitle(
                rf"Open Loop Frequency Response $\rightarrow$ {self.data_name}"
            )
            plt.savefig(
                Path(self.plots_dir)
                / f"open_loop_frequency_response_{self.data_name}.svg",
                format="svg",
            )

    def remove_delay(self, plot=True):
        omega = self.angular_frequencies
        f_obj_values = []
        tau_values = np.linspace(0, 0.25, 1000)
        for tau in tau_values:
            H_1 = self.H * np.exp(1j * omega * tau)
            sin_H1 = np.sin(np.angle(H_1))
            f_obj = np.sqrt(np.sum(sin_H1**2))
            f_obj_values.append(f_obj)

        f_obj_values = np.array(f_obj_values)

        if plot:
            plt.figure(figsize=(6, 5), dpi=180)
            plt.plot(tau_values, f_obj_values)
            plt.grid()
            plt.xlabel("Delay [s]")
            plt.ylabel("Objective Function Value")
            plt.tight_layout()

            plt.savefig(
                Path(self.plots_dir) / f"objetive_function_delay.svg",
                format="svg",
            )

        # Teste para tau > 0
        if plot:
            tau = 0.001
            H_1 = self.H * np.exp(1j * omega * tau)

            plt.figure(figsize=(6, 5), dpi=180)
            plt.semilogx(
                omega, np.rad2deg(np.unwrap(np.angle(self.H))), label="Original"
            )
            plt.semilogx(
                omega,
                np.rad2deg(np.unwrap(np.angle(H_1))),
                label=rf"Delay Removed (τ = {tau:.3g} s)",
            )
            plt.grid()
            plt.xlabel("Angular frequency [rad/s]")
            plt.ylabel("Phase [deg]")
            plt.tight_layout()
            plt.legend()

            plt.savefig(
                Path(self.plots_dir) / f"compare_original_delay_removed.svg",
                format="svg",
            )
