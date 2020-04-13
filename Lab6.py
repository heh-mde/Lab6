from random import *
from math import *
from numpy import *
from scipy.stats import f, t
from _pydecimal import Decimal
from datetime import datetime

# ConstLab5.py
p = 0.95  # Ймовірність для критеріїв
m = 3  # Початкова кількість експериментів
N = 4  # Початкова кількість дослідів
x1min = 15
x1max = 45
x2min = -25
x2max = 10
x3min = 45
x3max = 50
xlist = ((x1min, x1max), (x2min, x2max), (x3min, x3max))
rnd = 3  # Точність округлення
timing = 0 #Час виконання циклів

# Lab4.py
class lab5():
    def __init__(self, p, xlist, rnd):
        self.p = p
        self.q = round(1 - p, 2)  # Рівень значимості
        xcpmin = (xlist[0][0] + xlist[1][0] + xlist[2][0]) / 3
        xcpmax = (xlist[0][1] + xlist[1][1] + xlist[2][1]) / 3
        k = 3  # Кількість факторів
        self.l = sqrt(k)  # Значення плеча для k експериментів ПФЕ
        self.ymin = int(200 + xcpmin)
        self.ymax = int(200 + xcpmax)
        self.fisher = 0
        self.znach = []
        self.ydisplist = []
        self.xlist = xlist
        self.round = rnd

    def equation(self, N, m):
        rnd = self.round
        self.m = m
        self.ylist = []
        l = self.l
        xmatr = [[-1, -1, -1], [-1, -1, 1], [-1, 1, -1], [-1, 1, 1], [1, -1, -1], [1, -1, 1], [1, 1, -1], [1, 1, 1],
                 [-l, 0, 0], [l, 0, 0], [0, -l, 0], [0, l, 0], [0, 0, -l], [0, 0, l]]
        xmatr1 = xmatr
        x10 = (xlist[0][1] + xlist[0][0]) / 2
        x20 = (xlist[1][1] + xlist[1][0]) / 2
        x30 = (xlist[2][1] + xlist[2][0]) / 2
        x0 = [x10,x20,x30]
        for i in xmatr:
            k = -1
            for j in i:
                k += 1
                if j == 1:
                    i[k] = self.xlist[k][1]
                elif j == -1:
                    i[k] = self.xlist[k][0]
                else:
                    i[k] = i[k] * (self.xlist[k][1] - self.xlist[k][0]) / 2 + x0[k]
        for i in range(N):
            ylist1 = []
            for j in range(m):
                ylist1.append(round(
                    8.7 + 4.3 * xmatr[i][0] + 1.2 * xmatr[i][1] + 2.2 * xmatr[i][2] + 0.4 * xmatr[i][0] ** 2 + 1.0 *
                    xmatr[i][1] ** 2 + 6.4 * xmatr[i][2] ** 2 + 1.1 * xmatr[i][0] * xmatr[i][1] + 0.1 * xmatr[i][0] *
                    xmatr[i][2] + 9.2 * xmatr[i][1] * xmatr[i][2] + 1.2 * xmatr[i][0] * xmatr[i][1] * xmatr[i][
                        2] + randint(0,10) - 5, rnd))
            self.ylist.append(ylist1)
        ymed = []
        for i in range(N):
            ymed.append(round(sum(self.ylist[i]) / m, rnd))
        my = round(sum(ymed) / N, rnd)  # Коефіціент b0 норм.

        if N != 4:  # N = 8,14
            self.blist = []
            matrmlist = [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]]
            matrklist = [0, 0, 0, 0, 0, 0, 0, 0]
            if N == 14:
                matrklist.extend([0, 0, 0])
                matrmlist.extend([[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]])
                for i in range(11):
                    matrmlist[i].extend([0, 0, 0])
                for i in range(N):
                    xvarlist = (1, (xmatr1[i][0]), (xmatr1[i][1]), (xmatr1[i][2]), (xmatr1[i][0] * xmatr1[i][1]),
                                (xmatr1[i][0] * xmatr1[i][2]), (xmatr1[i][1] * xmatr1[i][2]),
                                (xmatr1[i][0] * xmatr1[i][1] * xmatr1[i][2]), xmatr1[i][0] * xmatr1[i][0],
                                xmatr1[i][1] * xmatr1[i][1], xmatr1[i][2] * xmatr1[i][2])
                    for j in range(len(matrmlist)):
                        for k in range(len(matrmlist)):
                            matrmlist[j][k] += xvarlist[j] * xvarlist[k]
                        matrklist[j] += ymed[i] * xvarlist[j]
                solve = linalg.solve(array(matrmlist), array(matrklist))  # Коефіціенти bi норм. N = 14
                for i in solve:
                    self.alist.append(round(i, rnd))

            for i in range(N):
                xvarlist = [1, (xmatr[i][0]), (xmatr[i][1]), (xmatr[i][2]), (xmatr[i][0] * xmatr[i][1]),
                            (xmatr[i][0] * xmatr[i][2]), (xmatr[i][1] * xmatr[i][2]),
                            (xmatr[i][0] * xmatr[i][1] * xmatr[i][2])]
                if N == 14:
                    xvarlist.extend([xmatr[i][0] * xmatr[i][0], xmatr[i][1] * xmatr[i][1], xmatr[i][2] * xmatr[i][2]])
                for j in range(len(matrmlist)):
                    for k in range(len(matrmlist)):
                        matrmlist[j][k] += xvarlist[j] * xvarlist[k]
                    matrklist[j] += ymed[i] * xvarlist[j]
            solve = linalg.solve(array(matrmlist), array(matrklist))  # Коефіціенти bi нaтур. N = 8,14

            if N == 8:
                a1 = round(  # Коефіціенти b1-b123 нормовані N = 8
                    (-1 * (ymed[0] + ymed[1] + ymed[2] + ymed[3]) + 1 * (ymed[4] + ymed[5] + ymed[6] + ymed[7])) / N,
                    rnd)
                a2 = round(
                    (-1 * (ymed[0] + ymed[1] + ymed[4] + ymed[5]) + 1 * (ymed[2] + ymed[3] + ymed[6] + ymed[7])) / N,
                    rnd)
                a3 = round(
                    (-1 * (ymed[0] + ymed[2] + ymed[4] + ymed[6]) + 1 * (ymed[1] + ymed[3] + ymed[5] + ymed[7])) / N,
                    rnd)
                a12 = round(
                    (-1 * (ymed[2] + ymed[3] + ymed[4] + ymed[5]) + 1 * (ymed[0] + ymed[1] + ymed[6] + ymed[7])) / N,
                    rnd)
                a13 = round(
                    (-1 * (ymed[1] + ymed[3] + ymed[4] + ymed[6]) + 1 * (ymed[0] + ymed[2] + ymed[5] + ymed[7])) / N,
                    rnd)
                a23 = round(
                    (-1 * (ymed[1] + ymed[2] + ymed[5] + ymed[6]) + 1 * (ymed[0] + ymed[1] + ymed[4] + ymed[7])) / N,
                    rnd)
                a123 = round(
                    (-1 * (ymed[0] + ymed[3] + ymed[5] + ymed[6]) + 1 * (ymed[1] + ymed[2] + ymed[4] + ymed[7])) / N,
                    rnd)
                self.alist = [my, a1, a2, a3, a12, a13, a23, a123]
            for i in solve:
                self.blist.append(round(i, rnd))
                self.xmatr = xmatr
        else:  # N = 4
            a1 = round((-1 * (ymed[0] + ymed[1]) + 1 * (ymed[2] + ymed[3])) / N, rnd)  # Коефіціенти b1-3 норм. N = 4
            a2 = round((-1 * (ymed[0] + ymed[2]) + 1 * (ymed[1] + ymed[3])) / N, rnd)
            a3 = round((-1 * (ymed[0] + ymed[3]) + 1 * (ymed[1] + ymed[2])) / N, rnd)

            deltax1 = (fabs(xlist[0][1] - xlist[0][0]) / 2)
            deltax2 = (fabs(xlist[1][1] - xlist[1][0]) / 2)
            deltax3 = (fabs(xlist[2][1] - xlist[2][0]) / 2)

            x10 = (xlist[0][1] + xlist[0][0]) / 2
            x20 = (xlist[1][1] + xlist[1][0]) / 2
            x30 = (xlist[2][1] + xlist[2][0]) / 2

            b0 = round(my - (a1 * x10 / deltax1) - (a2 * x20 / deltax2) - (a3 * x30 / deltax3),
                       2)  # Коефіціент b0 натур. N = 4
            b1 = round(a1 / deltax1, rnd)  # Коефіціенти b1-3 натур. N = 4
            b2 = round(a2 / deltax2, rnd)
            b3 = round(a3 / deltax3, rnd)
            self.blist = [b0, b1, b2, b3]
            self.alist = [my, a1, a2, a3]
            self.delta = [deltax1, deltax2, deltax3]
            self.x0 = [x10, x20, x30]
        self.ymed = ymed
        Task.printequa(N)

    def cochran(self, N, m):
        for i in range(N):
            ydisp = 0
            for k in range(m):
                ydisp += (self.ylist[i][k] - self.ymed[i]) ** 2
            ydisp /= m
            self.ydisplist.append(round(ydisp, self.round))

        self.groz = round(max(self.ydisplist) / sum(self.ydisplist), self.round)
        f1 = m - 1
        f2 = N
        partresult = self.q / f2
        params = [partresult, f1, (f2 - 1) * f1]
        fisher = f.isf(*params)
        result = fisher / (fisher + (f2 - 1))
        self.gkr = round(Decimal(result).quantize(Decimal('.0001')).__float__(), self.round)
        self.f1 = f1
        self.f2 = f2
        Task.printcoch(N)
        if self.groz < self.gkr:
            print("   Gp < Gкр => За критерієм Кохрана дисперсія однорідна з ймовірністю", p)
        else:
            print("   Gp > Gкр => За критерієм Кохрана дисперсія неоднорідна з ймовірністю", p)
            print("   Збільшуємо m на 1: m = {1}+1 = {0}".format(m + 1, m))
            m += 1
            Task.equation(N, m)
            Task.cochran(N, m)

    def student(self, N):
        self.d = N
        self.dvidtv = round(sum(self.ydisplist) / N, self.round)
        dkoef = round(self.d / (N * self.m), self.round)
        skoef = round(sqrt(dkoef), self.round)

        t0 = round(fabs(self.alist[0]) / skoef, self.round)
        t1 = round(fabs(self.alist[1]) / skoef, self.round)
        t2 = round(fabs(self.alist[2]) / skoef, self.round)
        t3 = round(fabs(self.alist[3]) / skoef, self.round)
        if N == 4:
            self.tlist = [t0, t1, t2, t3]
        else:
            t12 = round(fabs(self.alist[4]) / skoef, self.round)
            t13 = round(fabs(self.alist[5]) / skoef, self.round)
            t23 = round(fabs(self.alist[6]) / skoef, self.round)
            t123 = round(fabs(self.alist[7]) / skoef, self.round)
            self.tlist = [t0, t1, t2, t3, t12, t13, t23, t123]
            if N == 14:
                t11 = round(fabs(self.alist[8]) / skoef, self.round)
                t22 = round(fabs(self.alist[9]) / skoef, self.round)
                t33 = round(fabs(self.alist[10]) / skoef, self.round)
                self.tlist.extend([t11, t22, t33])

        self.f3 = self.f1 * self.f2
        tkr = Decimal(abs(t.ppf(self.q / 2, self.f3))).quantize(Decimal('.0001')).__float__()
        for troz in self.tlist:
            if troz < tkr:
                self.blist[self.tlist.index(troz)] = 0
                self.d -= 1
        blist = self.blist
        xlist = self.xlist
        if N != 4:
            self.yznachlist = []
            xmatr = self.xmatr
            for i in range(N):
                self.yznachlist.append(round(
                    blist[0] + blist[1] * xmatr[i][0] + blist[2] * xmatr[i][1] + blist[3] * xmatr[i][2] +
                    blist[4] * xmatr[i][0] * xmatr[i][1] + blist[5] * xmatr[i][0] * xmatr[i][2] + blist[6] *
                    xmatr[i][1] * xmatr[i][2] + blist[7] * xmatr[i][0] * xmatr[i][1] * xmatr[i][2], self.round))
            if N == 14:
                for i in range(N):
                    self.yznachlist[i] += round(
                        blist[8] * xmatr[i][0] ** 2 + blist[9] * xmatr[i][1] ** 2 + blist[10] * xmatr[i][2] ** 2,
                        self.round)
        else:
            y1 = round(
                blist[0] + blist[1] * xlist[0][0] + blist[2] * xlist[1][0] + blist[3] * xlist[2][0],
                self.round)
            y2 = round(
                blist[0] + blist[1] * xlist[0][0] + blist[2] * xlist[1][1] + blist[3] * xlist[2][1],
                self.round)
            y3 = round(
                blist[0] + blist[1] * xlist[0][1] + blist[2] * xlist[1][0] + blist[3] * xlist[2][1],
                self.round)
            y4 = round(
                blist[0] + blist[1] * xlist[0][1] + blist[2] * xlist[1][1] + blist[3] * xlist[2][0],
                self.round)
            self.yznachlist = [y1, y2, y3, y4]
        self.skoef = skoef
        self.dkoef = dkoef
        self.tkr = tkr
        Task.printstud(N)

    def fishers(self, N):
        dadekv = 0
        Nroz = N
        for k in range(N):
            dadekv += (self.yznachlist[k] - self.ymed[k]) ** 2
        if N == self.d:
            Nroz = N + 1
            self.Nroz = Nroz
        dadekv = round((dadekv * m) / (Nroz - self.d), self.round)
        froz = round(dadekv / self.dvidtv, self.round)
        self.f4 = Nroz - self.d
        fkr = Decimal(abs(f.isf(self.q, self.f4, self.f3))).quantize(Decimal('.0001')).__float__()
        self.dadekv = dadekv
        self.froz = froz
        self.fkr = fkr
        Task.printfish(N)
        if froz > fkr:
            if N == 14:
                print("   За критерієм Фішера рівняння з квадратичними членами неадекватне оригіналу з ймовірністю", p)
            else:
                print("   За критерієм Фішера лінійне рівняння неадекватне оригіналу з ймовірністю", p)
                if N == 4:
                    N = 8
                    print("\nПроведемо досліди для рівняння з ефектом взаємодії, тоді N =", N)
                elif N == 8:
                    N = 14
                    print("\nПроведемо досліди для рівняння з квадратичними членами, тоді N =", N)
                Task.equation(N, m)
                Task.cochran(N, m)
                Task.student(N)
                Task.fishers(N)
        else:
            print("   За критерієм Фішера лінійне рівняння регресії адекватне оригіналу з ймовірністю", p)

    def printequa(self, N):
        print("\nГенеруємо", N, "функцій відгуку для", self.m, "експериметнів:")
        for i in range(N):
            print("   Y{0} = {1}".format(i + 1, self.ylist[i]))
            aver = "Yсереднє" + str(i + 1) + " = ("
            for k in range(self.m - 1):
                aver += str(self.ylist[i][k]) + "+"
            aver += str(self.ylist[i][m - 1]) + ")/3 = " + str(self.ymed[i])
            print("  ", aver)
        if N == 4:
            print("\n1)Знайдемо коефіціенти рівняння регресії:\n\n  Для нормованих значень:")
            print("   b0 = my = ({0}+{1}+{2}+{3})/4 = {4}".format(*self.ymed, self.alist[0]))
            print("   b1 = a1 = (-1)*({0}+{1}) + 1*({2}+{3}) = {5}\n   b2 = a2 = (-1)*({0}+{2}) + 1*({1}+{3}) = {6}\n   \
b3 = a3 = (-1)*({0}+{3}) + 1*({1}+{2}) = {7}\n".format(*self.ymed, *self.alist))
            print("  Для натуральних:")
            for i in range(3):
                print("   Δx{3} = |({1})-({0})|/2 = {2}".format(*self.xlist[i], self.delta[i], i + 1))
                print("   x{3}0 = (({1})+({0}))/2 = {2}".format(*self.xlist[i], self.x0[i], i + 1))
            print("   b0 = {0} - ({1})*{4}/{7} - ({2})*{5}/{8} - ({3})*{6}/{6} - ({7})*{8}/{9} = {10}" \
                  .format(*self.alist, *self.x0, *self.delta, self.blist[0]))
            print("   b1 = {1}\n   b2 = {2}\n   b3 = {3}".format(*self.blist))
        elif N == 8:
            print("\n1)Знайдемо коефіціенти рівняння регресії:\n\n  Для нормованих значень:")
            print("   b0 = sum(ymed)\8 = {0}\n   b1 = sum(x1i*ymedi)\8 = {1}\n   b2  = sum(x2i*ymedi)\8 = {2}\n   "
                  "b3 = sum(x3i*ymedi)\8 = {3}\n   b12 = sum(x1i*x2i*ymedi)\8 = {4}\n   b13 = sum(x1i*x2i*ymedi)\8 = {5}"
                  "\n   b23 = sum(x1i*x2i*ymedi)\8 = {6}\n   b123 = sum(x1i*x2i*ymedi)\8 = {7}".format(*self.alist))
            print(" Для натуральних з системи рівнянь:")
            print("   b0 = {0}\n   b1 = {1}\n   b2 = {2}\n   b3 = {3}\n   b12 = {4}\n   b13 = {5}\n   b23 = {6}\n   "
                  "b123 = {7}   ".format(*self.blist))
        elif N == 14:
            print("\n1)Знайдемо коефіціенти рівняння регресії:\n\n  Для нормованих значень з системи рівнянь:")
            print("   b0 = {0}\n   b1 = {1}\n   b2 = {2}\n   b3 = {3}\n   b12 = {4}\n   b13 = {5}\n   b23 = {6}\n   "
                  "b123 = {7}\n   b11 = {8}\n   b22 = {9}\n   b33 = {10}".format(*self.alist))
            print(" \nДля натуральних з системи рівнянь:")
            print("   b0 = {0}\n   b1 = {1}\n   b2 = {2}\n   b3 = {3}\n   b12 = {4}\n   b13 = {5}\n   b23 = {6}\n   "
                  "b123 = {7}\n   b11 = {8}\n   b22 = {9}\n   b33 = {10}".format(*self.blist))

    def printcoch(self, N):
        print("\n2)Критерій Кохрана:\n\n  Знайдемо дисперсії по рядках:")
        for i in range(N):
            aver = "D" + "{y" + str(i + 1) + "} = ("
            for k in range(m - 1):
                aver += "(" + str(self.ylist[i][k]) + "-(" + str(self.ymed[i]) + "))^2 + "
            aver += "(" + str(self.ylist[i][m - 1]) + "-(" + str(self.ymed[i]) + "))^2)/" + str(m) + " = " + str(
                self.ydisplist[i])
            print("  ", aver)
        print(
            "\n   Dmax{{yi}} = {0}\n   Gp = {0}/({2}+{3}+{4}+...) = {1}".format(max(self.ydisplist),
                                                                                self.groz, *self.ydisplist))
        print("   f1 = {0} - 1 = {1}, f2 = 4, q = {3}\n   За таблицею Gкр = {2}".format(self.m, self.f1, self.gkr,
                                                                                        self.q))

    def printstud(self, N):
        if N == 4:
            print("\n2)Критерій Стьюдента:\n")
            print("   Dвідтворюваності = ({0}+{1}+{2}+{3})/4 = {4}".format(*self.ydisplist, self.dvidtv))
            print(
                "   D{{bi}} = {0}/(4*{1}) = {2}\n   S{{bi}} = sqrt({2}) = {3}".format(self.dvidtv, m, self.dkoef,
                                                                                      self.skoef))
            print("   t0 = |{0}|/{4} = {5}\n   t1 = |{1}|/{4} = {6}\n   t2 = |{2}|/{4} = {7}\n   \
t3 = |{3}|/{4} = {8}\n   ".format(*self.alist, self.skoef, *self.tlist))
            print("   f3 = 4*({0}-1) = {1}\n   За таблицею tkr = {2}".format(m, self.f3, self.tkr))
            for i in range(N):
                if self.tlist[i] < self.tkr:
                    print(
                        "   {0} < {1} => За критерієм Стьюдента коефіцієнт b{2} статистично незначущий з ймовірністю {3}"
                            .format(self.tlist[i], self.tkr, i, self.p))
                else:
                    print(
                        "   {0} > {1} => За критерієм Стьюдента коефіцієнт b{2} статистично значимий з ймовірністю {3}"
                            .format(self.tlist[i], self.tkr, i, self.p))
            print("\n   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*self.blist, x1min, x2min, x3min,
                                                                        self.yznachlist[0]))
            print("   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*self.blist, x1min, x2max, x3max,
                                                                      self.yznachlist[1]))
            print("   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*self.blist, x1max, x2min, x3max,
                                                                      self.yznachlist[2]))
            print("   {0} + {1}*{4} + {2}*{5} + {3}*{6} = {7}".format(*self.blist, x1max, x2max, x3min,
                                                                      self.yznachlist[3]))
        else:
            print("\n3)Критерій Стьюдента:\n")
            print("   Dвідтворюваності =", self.dvidtv)
            print(
                "   D{{bi}} =", self.dkoef, "\n   S{{bi}} =", self.skoef)
            for i in range(len(self.tlist)):
                print("   t{0} = {1}".format(i + 1, self.tlist[i]))
            print("   f3 = 4*({0}-1) = {1}\n   За таблицею tkr = {2}".format(m, self.f3, self.tkr))
            for i in range(len(self.tlist)):
                if self.tlist[i] < self.tkr:
                    print(
                        "   {0} < {1} => За критерієм Стьюдента коефіцієнт b{2} статистично незначущий з ймовірністю {3}"
                            .format(self.tlist[i], self.tkr, i, self.p))
                else:
                    print(
                        "   {0} > {1} => За критерієм Стьюдента коефіцієнт b{2} статистично значимий з ймовірністю {3}"
                            .format(self.tlist[i], self.tkr, i, self.p))
            print("")
            for i in range(N):
                print("   Y{0} = {1}".format(i + 1, self.yznachlist[i]))

    def printfish(self, N):
        print("\n3)Критерій Фішера:\n")
        print("   f4 = {2} - {0} = {1}".format(self.d, self.f4, self.Nroz))
        if N == 4:
            print(
                "   Dадекв =  {0}*(({5} - {1})**2 + ({6} - {2})**2 + ({7} - {3})**2 + ({8} - {4})**2)/({11}-{10}) = {9}"
                    .format(m, *self.ymed, *self.yznachlist, self.dadekv, self.d, self.Nroz))
        else:
            print("   Daдекв = ", self.dadekv)
        print("   Fр = {0}/{1} = {2}".format(self.dadekv, self.dvidtv, self.froz))
        print("   За таблицею Fкр =", self.fkr)


for i in range(100):
    s = datetime.now()
    Task = lab5(p, xlist, rnd)
    Task.equation(N, m)
    Task.cochran(N, m)
    Task.student(N)
    Task.fishers(N)
    timing += float(str(datetime.now() - s)[6:]) # Час виконання  100 циклів
print("\nСередній час виконання програми:",round(timing/100, 6),"секунди")
