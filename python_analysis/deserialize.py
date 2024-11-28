import struct
from dataclasses import dataclass

@dataclass
class Datapoint:
    t: float
    speed_x: float
    speed_y: float
    loc_x: float
    loy_y: float
    F_l: float

FIELD_COUNT = 6
T_EXPORTFLOAT_SIZE = 8
POINT_SIZE = FIELD_COUNT * T_EXPORTFLOAT_SIZE

def gen_fmt_string():
    str = "="
    for _ in range(FIELD_COUNT):
        str += "d"
    return str

FORMAT_STR = gen_fmt_string()
def read_point(idx, buffer):
    offset = idx * POINT_SIZE

    values = struct.unpack_from(FORMAT_STR, buffer, offset)
    return Datapoint(*values)

BUFFER_COUNT = 2 << 20
BUFFER_SIZE = BUFFER_COUNT * POINT_SIZE
class ExportReader:
    def __init__(self, path):
        self.file = open(path, "rb")
        self.buffer = None
        self.buffer_len = 0
        self.buffer_idx = 0
        self.fill_buffer()

    def destroy(self):
        self.file.close()
        self.file = None

    def fill_buffer(self):
        if not self.file.readable():
            self.buffer = None
            self.buffer_len = 0
            self.buffer_idx = 0
            return
        else:
            self.buffer = self.file.read(BUFFER_SIZE)
            buffer_len = len(self.buffer)
            if buffer_len == 0:
                self.buffer = None
                self.buffer_len = 0
                self.buffer_idx = 0
            else:
                self.buffer_len = len(self.buffer) / POINT_SIZE
        self.buffer_idx = 0

    def read_point(self):
        if self.buffer is None:
            return None
        point = read_point(self.buffer_idx, self.buffer)

        self.buffer_idx += 1
        if self.buffer_idx >= self.buffer_len:
            self.fill_buffer()

        return point