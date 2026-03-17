import matplotlib.pyplot as plt
import matplotlib.animation as animation


class Animator:
    def __init__(self, length, fps, updates_per_frame, ):
        self.length = length
        self.fps = fps
        self.updates_per_frame = updates_per_frame
