import math
import numpy as np
import matplotlib.pyplot as plt


# ------- Start T_Physics -------
class T_Physics:
    def __init__(self):
        self.g = 9.81
        self.air_density = 1.2

    def calc_gravity(self, a_mass):
        return self.g * a_mass * -1

    def calc_velocity(self, a_object, i):
        dt = a_object.t[i + 1] - a_object.t[i]
        dx = a_object.x[i + 1] - a_object.x[i]
        dy = a_object.y[i + 1] - a_object.y[i]
        v_x = dx / dt
        v_y = dy / dt
        return math.sqrt(v_x**2 + v_y**2)

    def calc_air_resistance(self, velocity):
        # preserving the sign of the force by using velocity * abs(velocity)
        return (
            -0.5
            * self.air_density
            * velocity
            * abs(velocity)
            * ball.C_d
            * ball.surface_area
        )

    def calc_energy(self, a_object):
        self.E_k = []
        self.E_p = []
        self.E_total = []
        for i in range(len(a_object.t) - 1):
            v = self.calc_velocity(a_object, i)
            self.E_k.append(0.5 * ball.mass * v**2)
            self.E_p.append(ball.mass * self.g * a_object.y[i])
            # calculate mechanical energy of last energy appended
            self.E_m = self.E_k[-1] + self.E_p[-1]
            self.E_total.append(self.E_m)


# ------- End of T_Physics -------


# ------- Start of T_Ball -------
class T_Ball:
    def __init__(self):
        self.radius = 0.0275
        self.mass = 0.06979
        self.surface_area = math.pi * self.radius**2
        self.C_d = 4.3


# ------- End of T_Ball -------


# ------- Start of T_Tracker_data -------
class T_Tracker_data:
    def __init__(self):
        self.data_values = open("tracker_data.txt", "r")
        self.data = self.data_values.read().split()
        self.data_values.close()
        self.t = []
        self.x = []
        self.y = []
        self.seperate_data()

    def get_v_0(self, value: str):
        selected_list = getattr(self, value)
        # using linear regression for the first 5 frames to estimate initial velocity more accurately
        t_values = self.t[:5]
        axis_values = selected_list[:5]
        # axis_values = a*t_values +b
        a, b = np.polyfit(t_values, axis_values, 1)
        return a

    def seperate_data(self):
        for i, value in enumerate(self.data):
            value = float(value)
            sort_logic = i % 3
            match sort_logic:
                case 0:
                    self.t.append(value)
                case 1:
                    self.x.append(value)
                case 2:
                    self.y.append(value)


# ------- End of T_Tracker_data -------


# ------- Start of T_math_model -------
class T_Math_model:
    def __init__(self, a_object):
        self.x = []
        self.y = []
        self.t = []
        self.F_g = physics.calc_gravity(a_object.mass)
        self.v_x = tracker_data.get_v_0("x")
        self.v_y = tracker_data.get_v_0("y")
        self.pos_x = tracker_data.x[0]
        self.pos_y = tracker_data.y[0]
        self.current_time = 0
        self.calculate_ball_movement(a_object)
        physics.calc_energy(self)

    def calculate_ball_movement(self, a_object):
        # loops untill ball hits "ground"
        while self.pos_y >= 0:
            # pushing all start values
            self.x.append(self.pos_x)
            self.y.append(self.pos_y)
            self.t.append(self.current_time)
            # chose small number for delta for more accurate values, 1000 positions pr second (1000 frames)
            dt = 0.001
            # calculate accelleration in each seperate direction
            a_x = physics.calc_air_resistance(self.v_x) / a_object.mass
            a_y = (physics.calc_air_resistance(self.v_y) + self.F_g) / a_object.mass
            self.current_time += dt
            self.v_x += a_x * dt
            self.v_y += a_y * dt
            self.pos_x += self.v_x * dt
            self.pos_y += self.v_y * dt


# ------- End of T_math_model -------


def draw_graphs(a_model):
    # --- Figure 1: Data -> Tracker ---
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(tracker_data.t, tracker_data.y, label="Tracker")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical position (m)")
    plt.title("Vertical position as a function of time")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(tracker_data.t, tracker_data.x)
    plt.xlabel("Time (s)")
    plt.ylabel("Horizontal position (m)")
    plt.title("Horizontal position as a function of time")
    plt.grid()
    plt.tight_layout()

    # --- Figure 2: Data -> Math model ---
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(a_model.t, a_model.y, label="Math model")
    plt.xlabel("Time (s)")
    plt.ylabel("Vertical position (m)")
    plt.title("Vertical position as a function of time")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(a_model.t, a_model.x, label="Math model")
    plt.xlabel("Time (s)")
    plt.ylabel("Horizontal position (m)")
    plt.title("Horizontal position as a function of time")
    plt.grid()
    plt.tight_layout()

    # --- Figure 3: Projectile Trajectory – Tracker vs model ---
    plt.figure()
    plt.plot(tracker_data.x, tracker_data.y, marker="o", label="Tracker")
    plt.plot(a_model.x, a_model.y, marker="x", label="Model")
    plt.title("Projectile Trajectory - Tracker vs model")
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.grid()
    plt.axis("equal")
    plt.legend()

    # --- Figure 4: Energy - Kinetic vs potential ---
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.plot(a_model.t[0 : len(physics.E_k)], physics.E_k)
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.title("Kinetic energy")
    plt.grid()

    plt.subplot(1, 2, 2)
    plt.plot(a_model.t[0 : len(physics.E_p)], physics.E_p)
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.title("Potential energy")
    plt.grid()
    plt.tight_layout()

    # --- Figure 5: Mechanical energy
    plt.figure()
    plt.plot(
        a_model.t[0 : len(physics.E_total)],
        physics.E_total,
        marker="o",
        label="Math model of ball",
    )
    plt.title("Mechanical energy")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.grid()
    plt.axis("equal")
    plt.legend()

    plt.show()


physics = T_Physics()
ball = T_Ball()
tracker_data = T_Tracker_data()
ball_model = T_Math_model(ball)
draw_graphs(ball_model)
