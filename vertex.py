class Vertex:
    def __init__(self, x, y, z=None, r=0, g=0, b=0):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.g = g
        self.b = b
        self.values = {}

    def attachValue(self, v_name, v_value):
        self.values[v_name] = v_value

    def getValue(self, v_name):
        return self.values.get(v_name)


class Vertex3D:
    def __init__(self, x, y, z, r, g, b):
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.g = g
        self.b = b
