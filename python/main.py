from generator import Generator as Gen
from generator import PriGenerateWay as Pgw
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.axes as axes

PRI_DETECTED_RANGE = [300, 700]  # pri detected range
TOA_RANGE = [0, 2*10**5]  # toa range

def sortingByTTP(toa: np.array,
                 bin_num: int) -> np.array:
    """
    find the pri for every toa

    """
    c = 0.9
    delta = 10

    def eraseNoiseByProb(mat: np.array) -> np.array:
        """
        Build histogram and delete noise

        """
        valid_pri = mat[np.nonzero(mat)]
        hist, bin_edges = np.histogram(valid_pri, bin_num, density=False)
        bin_width = (max_pri - min_pri) / bin_num
        mean_val = valid_pri.size / (max_pri - min_pri)
        for i, v in enumerate(hist):
            if v / bin_width <= c * mean_val:
                mat[np.logical_and(bin_edges[i] <= mat,  mat <= bin_edges[i + 1])] *= 0
        return mat

    min_pri = PRI_DETECTED_RANGE[0]
    max_pri = PRI_DETECTED_RANGE[1]

    # calculate toward direction and inverse direction plat transform matrix
    toward_ttp_mat = np.zeros((toa.size, toa.size), dtype=float)
    inverse_ttp_mat = np.zeros((toa.size, toa.size), dtype=float)

    for k in range(toa.size):
        index = 0
        for j in range(k + 1, toa.size, 1):
            cur_pri = toa[j] - toa[k]
            if min_pri <= cur_pri <= max_pri:
                toward_ttp_mat[k][index] = cur_pri
                index += 1

    for k in range(toa.size):
        index = 0
        for j in range(k - 1, -1, -1):
            cur_pri = toa[k] - toa[j]
            if min_pri <= cur_pri <= max_pri:
                inverse_ttp_mat[k][index] = cur_pri
                index += 1

    # sort pri by probability and delete noise
    toward_ttp_mat = eraseNoiseByProb(toward_ttp_mat)
    inverse_ttp_mat = eraseNoiseByProb(inverse_ttp_mat)

    # relative analysis
    for k in range(toa.size):
        for i in range(toa.size):
            u = toward_ttp_mat[k][i]
            v = inverse_ttp_mat[k][i]
            if u and v:
                if abs(v - u) > delta:
                    toward_ttp_mat[k][i] = 0
            else:
                break

    # delete mirror points
    ttp_mat = np.zeros(toa.size)
    for k in range(toa.size):
        row = toward_ttp_mat[k, :]
        val = row[np.nonzero(row)]
        if val.size:
            ttp_mat[k] = val[0]

    return ttp_mat

def plotTTPMat(ax: axes,
               toa: np.array,
               ttp: np.array,
               title: str = None,
               plot_ylabel: bool = False) -> None:
    ax.scatter(toa, ttp, marker='o', c=ttp, alpha=0.75)
    ax.set_ylim(PRI_DETECTED_RANGE[0], PRI_DETECTED_RANGE[1])
    ax.xaxis.get_major_formatter().set_powerlimits((0, 1))
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('toa[us]', fontsize=8, loc='right')
    if plot_ylabel:
        ax.set_ylabel('pri[us]', fontsize=9)

def sortingByPRITransform(toa: np.array,
                          bin_num: int) -> tuple:
    bins = np.linspace(PRI_DETECTED_RANGE[0], PRI_DETECTED_RANGE[1], num=bin_num + 1, endpoint=True)
    d = np.zeros(bin_num, dtype=complex)
    for k in range(len(bins) - 1):
        for i in range(len(toa)):
            for j in range(i + 1, len(toa), 1):
                pri = toa[j] - toa[i]
                if bins[k] <= pri <= bins[k + 1]:
                    d[k] = d[k] + np.exp(2 * np.pi * toa[j] * 1j / (toa[j] - toa[i]))
                elif pri > bins[k + 1]:
                    break

    bins += (PRI_DETECTED_RANGE[1] - PRI_DETECTED_RANGE[0]) / bin_num / 2
    return bins[:-1], abs(d)

def plotPRITransformVal(ax: axes,
                        bins: np.array,
                        vals: np.array,
                        title: str = None,
                        plot_ylabel: bool = False) -> None:
    ax.plot(bins, vals, linewidth=2)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('pri[us]', fontsize=8, loc='right')
    if plot_ylabel:
        ax.set_ylabel('interval values', fontsize=10)

def run() -> None:
    # generate toa by five ways
    stable_toa1 = Gen.generateTOA([400], TOA_RANGE, Pgw.STABLE)
    stable_toa2 = Gen.generateTOA([610], TOA_RANGE, Pgw.STABLE)
    stable_toa = mergeSortedArray(stable_toa1, stable_toa2)

    wobble_toa1 = Gen.generateTOA([400, 0.3], TOA_RANGE, Pgw.WOBBLE)
    wobble_toa2 = Gen.generateTOA([610, 0.3], TOA_RANGE, Pgw.WOBBLE)
    wobble_toa = mergeSortedArray(wobble_toa1, wobble_toa2)

    irregular_toa1 = Gen.generateTOA([400, 430, 450], TOA_RANGE, Pgw.IRREGULAR)
    irregular_toa2 = Gen.generateTOA([620, 640, 650], TOA_RANGE, Pgw.IRREGULAR)
    irregular_toa = mergeSortedArray(irregular_toa1, irregular_toa2)

    slip_toa1 = Gen.generateTOA([400, 450, 5], TOA_RANGE, Pgw.SLIP)
    slip_toa2 = Gen.generateTOA([610, 670, 10], TOA_RANGE, Pgw.SLIP)
    slip_toa = mergeSortedArray(slip_toa1, slip_toa2)

    sin_toa1 = Gen.generateTOA([400, 0.05, 50, 0], TOA_RANGE, Pgw.SIN)
    sin_toa2 = Gen.generateTOA([610, 0.07, 70, 0], TOA_RANGE, Pgw.SIN)
    sin_toa = mergeSortedArray(sin_toa1, sin_toa2)

    '''
    # estimate pri
    stable_pri = sortingByTTP(stable_toa, bin_num=10)
    wobble_ttp = sortingByTTP(wobble_toa, bin_num=10)
    irregular_pri = sortingByTTP(irregular_toa, bin_num=10)
    slip_pri = sortingByTTP(slip_toa, bin_num=10)
    sin_pri = sortingByTTP(sin_toa, bin_num=10)

    # plot toa - ttp distribution
    fig, axs = plt.subplots(2, 3, sharey='row')
    plotTTPMat(axs[0, 0], stable_toa, stable_pri, plot_ylabel=True, title='Two Stable Pri TTP Transform')

    plotTTPMat(axs[0, 1], wobble_toa, wobble_ttp, title='Two Wobble Pri TTP Transform')

    plotTTPMat(axs[0, 2], irregular_toa, irregular_pri, title='Two irregular Pri Transform')

    plotTTPMat(axs[1, 0], slip_toa, slip_pri, title='Two Slip Pri Transform')

    plotTTPMat(axs[1, 1], sin_toa, sin_pri, plot_ylabel=True, title='Two Close Sin Pri Transform')

    plt.show()
    '''

    # estimate pri
    stable_pri, stable_d = sortingByPRITransform(stable_toa, bin_num=200)
    wobble_pri, wobble_d = sortingByPRITransform(wobble_toa, bin_num=200)
    irregular_pri, irregular_d = sortingByPRITransform(irregular_toa, bin_num=200)
    slip_pri, slip_d = sortingByPRITransform(slip_toa, bin_num=200)
    sin_pri, sin_d = sortingByPRITransform(sin_toa, bin_num=200)

    # plot toa - ttp distribution
    fig, axs = plt.subplots(2, 3, sharey='row')
    plotPRITransformVal(axs[0, 0], stable_pri, stable_d, plot_ylabel=True, title='Two Stable Pri TTP Transform')

    plotPRITransformVal(axs[0, 1], wobble_pri, wobble_d, title='Two Wobble Pri TTP Transform')

    plotPRITransformVal(axs[0, 2], irregular_pri, irregular_d, title='Two irregular Pri Transform')

    plotPRITransformVal(axs[1, 0], slip_pri, slip_d, plot_ylabel=True, title='Two Slip Pri Transform')

    plotPRITransformVal(axs[1, 1], sin_pri, sin_d, title='Two Close Sin Pri Transform')

    plt.show()


def mergeSortedArray(arr1: np.array,
                     arr2: np.array) -> np.array:
    size = len(arr1) + len(arr2)
    ret_arr = np.empty(size)
    k1 = 0; k2 = 0; k = 0
    while k1 < len(arr1) or k2 < len(arr2):
        if k1 < len(arr1) and k2 < len(arr2):
            if arr1[k1] < arr2[k2]:
                ret_arr[k] = arr1[k1]
                k1 += 1
            else:
                ret_arr[k] = arr2[k2]
                k2 += 1
        elif k1 < len(arr1):
            ret_arr[k] = arr1[k1]
            k1 += 1
        else:
            ret_arr[k] = arr2[k2]
            k2 += 1
        k += 1
    return ret_arr

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    run()
