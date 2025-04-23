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
    def __init__(self, a_object, air_resistance):
        self.x = []
        self.y = []
        self.t = []
        self.F_g = physics.calc_gravity(a_object.mass)
        self.v_x = tracker_data.get_v_0("x")
        self.v_y = tracker_data.get_v_0("y")
        self.current_time = 0
        self.calculate_ball_movement(a_object, air_resistance)
        """
        if air_resistance:
            physics.calc_energy(self)
        """

    def calculate_ball_movement(self, a_object, air_resistance):
        self.pos_x = tracker_data.x[0]
        self.pos_y = tracker_data.y[0]
        # pushing all start values
        self.x.append(self.pos_x)
        self.y.append(self.pos_y)
        self.t.append(self.current_time)
        # chose small number for delta for more accurate values, 1000 positions pr second (1000 frames)
        dt = 0.001
        # loops untill ball hits "ground"
        while self.pos_y >= 0:
            self.current_time += dt
            if air_resistance:
                # calculate accelleration in each seperate direction
                a_x = physics.calc_air_resistance(self.v_x) / a_object.mass
                a_y = (physics.calc_air_resistance(self.v_y) + self.F_g) / a_object.mass
            else:
                print("heisann")
                a_x = 0
                a_y = physics.g * -1

            self.v_x += a_x * dt
            self.v_y += a_y * dt
            self.pos_x += self.v_x * dt
            self.pos_y += self.v_y * dt

            self.x.append(self.pos_x)
            self.y.append(self.pos_y)
            self.t.append(self.current_time)


# ------- End of T_math_model -------


def draw_graphs(model_with_Fd, model_without_Fd):
    # --- Figure 1: Projectile Trajectory – Tracker vs model ---
    plt.figure()
    plt.plot(tracker_data.x, tracker_data.y, marker="o", label="Measurements")
    plt.plot(model_with_Fd.x, model_with_Fd.y, marker=",", label="Model")
    plt.plot(
        model_without_Fd.x, model_without_Fd.y, marker=",", label="Model without Fd"
    )
    plt.title("Projectile Trajectory - Measurements vs model")
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.grid()
    plt.axis("equal")
    plt.legend()

    # --- Figure 2: Kinetic energy - with Fd  vs wihtout Fd ---
    plt.figure()
    physics.calc_energy(model_with_Fd)
    plt.plot(
        model_with_Fd.t[0 : len(physics.E_k)],
        physics.E_k,
        marker=",",
        label="Model",
    )

    physics.calc_energy(model_without_Fd)
    plt.plot(
        model_without_Fd.t[0 : len(physics.E_k)],
        physics.E_k,
        marker=",",
        label="Model without Fd",
    )
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.title("Kinetic energy")
    plt.grid()
    plt.axis("equal")
    plt.legend()

    # --- Figure 3: Potential energy - with Fd vs without Fd ---
    plt.figure()
    physics.calc_energy(model_with_Fd)
    plt.plot(
        model_with_Fd.t[0 : len(physics.E_p)],
        physics.E_p,
        marker=",",
        label="Model",
    )
    physics.calc_energy(model_without_Fd)
    plt.plot(
        model_without_Fd.t[0 : len(physics.E_p)],
        physics.E_p,
        marker=",",
        label="Model without Fd",
    )
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.title("Potential energy")
    plt.grid()
    plt.axis("equal")
    plt.legend()
    # --- Figure 3: Mechanical energy
    plt.figure()
    physics.calc_energy(model_with_Fd)
    plt.plot(
        model_with_Fd.t[0 : len(physics.E_total)],
        physics.E_total,
        marker=",",
        label="Model",
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
ball_model = T_Math_model(ball, True)
ball_model_no_Fd = T_Math_model(ball, False)
draw_graphs(ball_model, ball_model_no_Fd)
