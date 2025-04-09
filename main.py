import math
import numpy as np
import matplotlib.pyplot as plt


data_values = open("tracker_data.txt", "r")
data = data_values.read().split()


# ------- Start of T_Tracker_data -------
class T_Tracker_data:
    def __init__(self):
        self.t = []
        self.x = []
        self.y = []
        self.positions = []
        self.seperate_data()

    def get_v_0(self, value: str):
        selected_list = getattr(self, value)
        """
        dt = self.t[1] - self.t[0]
        return (selected_list[1] - selected_list[0]) / dt
        """

        # using linear regression for the first 5 frames to estimate initial velocity more accurately
        t_values = self.t[:5]
        axis_values = selected_list[:5]
        # axis_values = a*t_values +b
        a, b = np.polyfit(t_values, axis_values, 1)
        return a

    def seperate_data(self):
        for i, value in enumerate(data):
            value = float(value)
            sort_logic = i % 3
            match sort_logic:
                case 0:
                    self.t.append(value)
                case 1:
                    self.x.append(value)
                case 2:
                    self.y.append(value)

    def estimate_drag_coefficient():
        # asdf

        return C_d


# ------- End of T_Tracker_data -------


# ------- Start of T_math_model -------
class T_math_model:
    def __init__(self):
        self.x = []
        self.y = []
        self.t = []
        # mass in kg (69.79g)
        self.m = 0.06979
        # radius in meter (5.5 cm)
        self.r = 0.0275
        self.air_density = 1.225
        self.F_g = (9.81 * self.m) * -1
        # C_d is a temp value, remember to calcualte the real drag-coefficient
        self.C_d = 2.3
        # self.C_d = tracker_data.estimate_drag_coefficient()
        self.v_x = tracker_data.get_v_0("x")
        self.v_y = tracker_data.get_v_0("y")
        self.pos_x = tracker_data.x[0]
        self.pos_y = tracker_data.y[0]
        self.current_time = 0
        self.calculate_ball_movement()

    # fmt: off
    def calc_air_resistance(self, velocity):
        # preserving the sign of the force by using velocity * abs(velocity)
        return (-0.5*self.air_density*velocity*abs(velocity)*self.C_d*math.pi*self.r**2)
    # fmt: on

    def calculate_ball_movement(self):
        # loops untill ball hits "ground"
        while self.pos_y >= 0:
            # pushing all start values
            self.x.append(self.pos_x)
            self.y.append(self.pos_y)
            self.t.append(self.current_time)
            # chose small number for delta for more accurate values, 1000 positions pr second (1000 frames)
            dt = 0.001
            # calculate accelleration in each seperate direction
            a_x = self.calc_air_resistance(self.v_x) / self.m
            a_y = (self.calc_air_resistance(self.v_y) + self.F_g) / self.m
            self.current_time += dt
            self.v_x += a_x * dt
            self.v_y += a_y * dt
            self.pos_x += self.v_x * dt
            self.pos_y += self.v_y * dt

    def calc_energy(self):
        E_k = None
        E_p = None
        for i in range(len(self.t)):
            v = math.sqrt((self.x[i]) ** 2 + (self.y[i]) ** 2)
            E_k = 0.5 * self.m * v**2
            E_p = self.m * 9.91 * self.y[i]
        E_m = E_k + E_p
        print(E_m)


# ------- End of T_math_model -------


def draw_graphs():
    # --- Figure 1: Data -> Tracker ---
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(tracker_data.t, tracker_data.y, label="Tracker")
    plt.xlabel("time (s)")
    plt.ylabel("y - position (m)")
    plt.title("y as a function of time")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(tracker_data.t, tracker_data.x)
    plt.xlabel("time (s)")
    plt.ylabel("x - position (m)")
    plt.title("x as a function of time")
    plt.grid()
    plt.tight_layout()

    # --- Figure 2: Data -> Math model ---
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(test.t, test.y, label="Math model")
    plt.xlabel("time (s)")
    plt.ylabel("y - position (m)")
    plt.title("y as a function of time")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(test.t, test.x, label="Math model")
    plt.xlabel("time (s)")
    plt.ylabel("x - position (m)")
    plt.title("x as a function of time")
    plt.grid()
    plt.tight_layout()

    # === Figure 5: Kastebane – Tracker vs Modell ===
    plt.figure()
    plt.plot(tracker_data.x, tracker_data.y, marker="o", label="Tracker")
    plt.plot(test.x, test.y, marker="x", label="Modell")
    plt.title("Kastebane – Tracker vs Modell")
    plt.xlabel("x - position (m)")
    plt.ylabel("y - position (m)")
    plt.grid()
    plt.axis("equal")
    plt.legend()

    plt.show()


tracker_data = T_Tracker_data()
test = T_math_model()
test.calc_energy()
draw_graphs()
# ---- TO DO ----
#   * beregne C_d og bytte ut med 0.47
#   * plot (tegn graf) for energi over tid
