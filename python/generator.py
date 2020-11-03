import math
import numpy as np
from enum import Enum
import random

class PriGenerateWay(Enum):
    STABLE = 1
    WOBBLE = 2
    IRREGULAR = 3
    SLIP = 4
    SIN = 5

class Generator:
    toa = []

    @staticmethod
    def generateTOA(pri_params: list[float],
                    toa_range: list[float],
                    pri_generate_way: PriGenerateWay,
                    float_rate: float = 0.01,
                    reject_rate: float = 0.1) -> np.array:
        """
        Generate toa.

        """
        pri_generator = None
        if pri_generate_way == PriGenerateWay.STABLE:
            pri_generator = Generator.generateStablePRI
        elif pri_generate_way == PriGenerateWay.WOBBLE:
            pri_generator = Generator.generateWobblePRI
        elif pri_generate_way == PriGenerateWay.IRREGULAR:
            pri_generator = Generator.generateIrregularPRI
        elif pri_generate_way == PriGenerateWay.SLIP:
            pri_generator = Generator.generateSlipPRI
        elif pri_generate_way == PriGenerateWay.SIN:
            pri_generator = Generator.generateSinPRI

        Generator.toa = []
        toa_begin = toa_range[0]
        toa_end = toa_range[1]
        Generator.toa.append(toa_begin)
        rng = np.random.default_rng()

        while Generator.toa[-1] <= toa_end:
            current_pri = pri_generator(pri_params)
            Generator.toa.append(Generator.toa[-1] + (1 + (2 * rng.random() - 1) * float_rate) * current_pri)

        Generator.toa.pop()
        reserved_list = random.sample(range(len(Generator.toa)), math.floor((1 - reject_rate) * len(Generator.toa)))
        toa_arr = np.fromiter(Generator.toa, float)
        reserved_list.sort()
        return toa_arr[reserved_list]

    @staticmethod
    def generateStablePRI(params: list[float]) -> float:
        """
        Generate stable pri.

        """
        return params[0]

    @staticmethod
    def generateWobblePRI(params: list[float]) -> float:
        """
        Generate wobble pri.
        """
        pri = params[0]
        wobble_rate = params[1]
        rng = np.random.default_rng()
        return pri * (1 + math.sqrt(1 / 3) * rng.standard_normal() * wobble_rate)

    @staticmethod
    def generateIrregularPRI(params: list[float]) -> float:
        """
        Generate wobble pri.

        """
        k = (len(Generator.toa) - 1) % len(params)
        pri = params[k]
        return pri

    @staticmethod
    def generateSlipPRI(params: list[float]) -> float:
        """
        Generate slip pri.

        """
        if len(Generator.toa) < 2:
            return params[0]
        else:
            pri = Generator.toa[-1] - Generator.toa[-2] + params[2]
            if pri > params[1]:
                return params[0]
            else:
                return pri

    @staticmethod
    def generateSinPRI(params: list[float]) -> float:
        """
        Generate sin pri.

        """
        mid_val = params[0]
        am = params[1]
        interval = params[2]
        ph = params[3]
        return mid_val + am * mid_val * np.sin(2 * np.pi * (len(Generator.toa) / interval) + ph)




