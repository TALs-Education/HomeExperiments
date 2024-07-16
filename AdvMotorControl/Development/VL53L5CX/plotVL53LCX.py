import serial
import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

class SerialReader:
    def __init__(self, port, baud_rate=9600, grid_size=(4, 4)):
        self.port = port
        self.baud_rate = baud_rate
        self.grid_size = grid_size
        self.ser = None

    def connect(self):
        try:
            self.ser = serial.Serial(self.port, self.baud_rate)
            print(f"Connected to {self.port} at {self.baud_rate} baud rate.")
        except serial.SerialException as e:
            print(f"Failed to connect to {self.port}: {e}")
            return False
        return True

    def read_line(self):
        if self.ser and self.ser.in_waiting > 0:
            return self.ser.readline().decode('utf-8').rstrip()
        return None

    def close(self):
        if self.ser:
            self.ser.close()
            print(f"Disconnected from {self.port}")

    def start_plot(self):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlabel('X')
        self.ax.set_ylabel('Y')
        self.ax.set_zlabel('Z')

        # Set up the animation
        self.anim = FuncAnimation(self.fig, self.update_plot, interval=50, blit=False)
        plt.show()

    def update_plot(self, frame):
        if not self.ser.in_waiting:
            return

        line = self.read_line()
        print(line)  # Echo the data to the terminal
        
        if line:
            values = [int(value) for value in line.split(',') if value.strip()]
            if len(values) == self.grid_size[0] * self.grid_size[1]:
                Z = np.array(values).reshape(self.grid_size)
                x = np.linspace(0, self.grid_size[0]-1, self.grid_size[0])
                y = np.linspace(0, self.grid_size[1]-1, self.grid_size[1])
                X, Y = np.meshgrid(x, y)

                self.ax.clear()
                self.ax.set_xlabel('X')
                self.ax.set_ylabel('Y')
                self.ax.set_zlabel('Z')
                self.ax.plot_surface(X, Y, Z.T, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
                self.ax.set_zlim(0, 4000)  # Set the z-axis range

# Example usage
if __name__ == "__main__":
    serial_port = 'COM3'
    baud_rate = 9600
    grid_size = (4, 4)

    reader = SerialReader(serial_port, baud_rate, grid_size)

    if reader.connect():
        reader.start_plot()
        reader.close()
