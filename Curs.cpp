#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <cmath>
using namespace std;

// Вспомогательная функция для сравнения float
bool floatEqual(float a, float b) {
    return  abs(a - b) < 1e-6f;
}

// Функция для сравнения двух точек
bool pointEqual(float* p1, float* p2) {
    return floatEqual(p1[0], p2[0]) && floatEqual(p1[1], p2[1]);
}

// Функция для копирования точки
void copyPoint(float* dest, float* src) {
    dest[0] = src[0];
    dest[1] = src[1];
}

// Функция для сравнения двух шестиугольников (лексикографически)
int compareHexagons(float* h1[6], float* h2[6]) {
    for (int i = 0; i < 6; i++) {
        if (!floatEqual(h1[i][0], h2[i][0])) {
            return (h1[i][0] < h2[i][0]) ? -1 : 1;
        }
        if (!floatEqual(h1[i][1], h2[i][1])) {
            return (h1[i][1] < h2[i][1]) ? -1 : 1;
        }
    }
    return 0; // равны
}

// Функция для получения канонической формы шестиугольника
void getCanonicalForm(float* hex[6], float* result[6]) {
    float* permutations[6][6]; // 6 перестановок по 6 указателей на точки

    // Генерируем все циклические перестановки
    for (int shift = 0; shift < 6; shift++) {
        for (int i = 0; i < 6; i++) {
            int idx = (i + shift) % 6;
            permutations[shift][i] = hex[idx];
        }
    }

    // Находим минимальную перестановку
    int minIdx = 0;
    for (int shift = 1; shift < 6; shift++) {
        if (compareHexagons(permutations[shift], permutations[minIdx]) < 0) {
            minIdx = shift;
        }
    }

    // Копируем минимальную перестановку в результат
    for (int i = 0; i < 6; i++) {
        result[i] = permutations[minIdx][i];
    }
}

// Функция обмена двух шестиугольников в массиве
void swapHexagons(float*** points, int i, int j) {
    float** temp = points[i];
    points[i] = points[j];
    points[j] = temp;
}

// Функция для получения канонической формы шестиугольника по индексу
void getCanonicalFormAtIndex(float*** points, int index, float* result[6]) {
    getCanonicalForm(points[index], result);
}

// Упрощенная версия сортировки (без выделения дополнительной памяти)
void sortHexagons(float*** points, int size) {
    // Пузырьковая сортировка на основе канонических форм
    for (int i = 0; i < size - 1; i++) {
        for (int j = 0; j < size - i - 1; j++) {
            float* canonical1[6];
            float* canonical2[6];

            getCanonicalForm(points[j], canonical1);
            getCanonicalForm(points[j + 1], canonical2);

            if (compareHexagons(canonical1, canonical2) > 0) {
                swapHexagons(points, j, j + 1);
            }
        }
    }
}

// Функция для группировки циклических перестановок
void groupCyclicPermutations(float*** points, int size) {
    sortHexagons(points, size);
}

// Функция для вывода шестиугольника (для отладки)
void printHexagon(float* hex[6]) {
    for (int i = 0; i < 6; i++) {
         cout << "(" << hex[i][0] << ", " << hex[i][1] << ") ";
    }
     cout <<  endl;
}

// Функция для проверки, является ли один шестиугольник циклической перестановкой другого
bool isCyclicPermutation(float* h1[6], float* h2[6]) {
    // Проверяем все возможные сдвиги
    for (int shift = 0; shift < 6; shift++) {
        bool match = true;
        for (int i = 0; i < 6; i++) {
            int idx = (i + shift) % 6;
            if (!pointEqual(h1[i], h2[idx])) {
                match = false;
                break;
            }
        }
        if (match) return true;
    }
    return false;
}


void srt(float*** she, float*** ps) {
    //вот тут будет функция, чтобы создать оригинальные шестерки
}

float mmm(float* M, int  r) {
    int p = -1;
    for (int i = 0; i < r; i++) {
        if (M[i] > p) {
            p = M[i];
        }
    }
    return p;
}

void pr_m(float*** M, ofstream& log, int c, ofstream& out, float* pro, int q) {
    //M - массив шестерок
    //с - наибольшее кол-во элементов в шестиугольнике
    log << "\nНаилучшие шестерки: \n";
    out << "\nНаилучшие шестерки: \n";
    for (int p = 0; p < q; p++) {
        if (pro[p] == c) {
            for (int i = 0; i < 6; i++) {
                float a = M[p][i][0];
                float z = M[p][i][1];
                log << a << " " << z << "\n";
                out << a << " " << z << "\n";
            }
            log << "-----\n";
            out << "-----\n";
        }

    }
    log << "Внутри шестиугольников " << c << " элементов";
    out << "Внутри шестиугольников " << c << " элементов";
}

bool prin(float*** she, int j, float ax, float ay, float s, float p) {
    float ss = 0;
    for (int i = 0; i < 5; i++) {
        float ax1 = she[j][i][0];
        float ay1 = she[j][i][1];
        float ax2 = she[j][i + 1][0];
        float ay2 = she[j][i + 1][1];

        ss += ((abs((ax2 - ax1) * (ay - ay1) - (ax - ax1) * (ay2 - ay1))) / 2);
    }
    float aax1 = she[j][0][0];
    float aay1 = she[j][0][1];
    float aax2 = she[j][5][0];
    float aay2 = she[j][5][1];

    ss += ((abs(((aax2 - aax1) * (ay - aay1)) - ((ax - aax1) * (aay2 - aay1)))) / 2);

    if ((abs(ss - s)) <= p * s) return 1;
    return 0;
}

float pl(float*** p, int n) {
    float ax1 = p[n][0][0];
    float ay1 = p[n][0][1];
    float ax2 = p[n][1][0];
    float ay2 = p[n][1][1];
    float a = sqrt((pow((ax1 - ax2), 2)) + pow((ay1 - ay2), 2));
    float S = ((3 * sqrt(3) * a * a) / 2);
    return S;
}

void kal(float** p, float* p_pro, float*** she, int ss, int m, float pp, ofstream& log) {
    log << "\nЗаход в kal\n";

    //начинаем проходиться по всем точкам
    for (int i = 0; i < m; i++) {
        float ax1 = p[i][0];
        float ay1 = p[i][1];
        //заходим в одну из шестерок
        for (int j = 0; j < ss; j++) {
            bool flag = true;

            //проходимся по эл-там шестерки
            for (int jj = 0; jj < 6; jj++) {
                float ax2 = she[j][jj][0];
                float ay2 = she[j][jj][1];

                if ((abs(ax1 - ax2) < pp) && (abs(ay1 - ay2) < pp)) flag = false;
            }
            if (flag) {
                //окей, точка не является частью шестерки, так что проверяем ее принадлежность
                float S = pl(she, j);

                if (prin(she, j, ax1, ay1, S, pp)) {

                    log << "\nN: " << j << " |------------подходит: " << ax1 << " " << ay1;
                    p_pro[j]++;
                }
                else log << "\nN: " << j << " |НЕ подходит: " << ax1 << " " << ay1;
            }
        }
    }
}

bool v_t1_t2(float M[6][2], float p) {

    bool a1 = ((M[2][0] - M[1][0]) * (M[3][1] - M[2][1]) - (M[3][0] - M[2][0]) * (M[2][1] - M[1][1])) > 0;
    //(x2-x1)(y3-y2) - (x3-x2)(y2-y1) > 0
    bool a2 = ((M[3][0] - M[2][0]) * (M[4][1] - M[3][1]) - (M[4][0] - M[3][0]) * (M[3][1] - M[2][1])) > 0;
    //(x3-x2)(y4-y3) - (x4-x3)(y3-y2) > 0
    bool a3 = ((M[1][0] - M[0][0]) * (M[2][1] - M[1][1]) - (M[2][0] - M[1][0]) * (M[1][1] - M[0][1])) > 0;
    //(x1 - x0)(y2 - y1) - (x2 - x1)(y1 - y0) > 0

    if (a1 && a2 && a3) return true;
    else return 0;
}

bool ds_t1_t2(float M[6][2], float p) {
    float a1 = sqrt(pow((M[0][0] - M[1][0]), 2) + pow((M[0][1] - M[1][1]), 2));
    float a2 = sqrt(pow((M[1][0] - M[2][0]), 2) + pow((M[1][1] - M[2][1]), 2));
    float a3 = sqrt(pow((M[2][0] - M[3][0]), 2) + pow((M[2][1] - M[3][1]), 2));
    float a4 = sqrt(pow((M[3][0] - M[4][0]), 2) + pow((M[3][1] - M[4][1]), 2));
    float a5 = sqrt(pow((M[4][0] - M[5][0]), 2) + pow((M[4][1] - M[5][1]), 2));
    float a6 = sqrt(pow((M[5][0] - M[0][0]), 2) + pow((M[5][1] - M[0][1]), 2));
    if (a1 == a2 && a1 == a3 && a1 == a4 && a1 == a5 && a1 == a6) return true;
    else if (((sqrt(pow((a1 - a2), 2))) < p) && ((sqrt(pow((a1 - a3), 2))) < p) && ((sqrt(pow((a1 - a4), 2))) < p) && ((sqrt(pow((a1 - a5), 2))) < p) && ((sqrt(pow((a1 - a6), 2))) < p)) return true;
    else return false;
}

bool c_t1_t2(float M[6][2], float p) {
    float ax1 = (M[0][0] + M[3][0]) / 2;
    float ay1 = (M[0][1] + M[3][1]) / 2;

    float ax2 = (M[2][0] + M[5][0]) / 2;
    float ay2 = (M[2][1] + M[5][1]) / 2;

    float ax3 = (M[1][0] + M[4][0]) / 2;
    float ay3 = (M[1][1] + M[4][1]) / 2;
    if ((ax1 == ax2 && ax1 == ax3) && (ay1 == ay2 && ay1 == ay3)) return true;
    else if ((((sqrt(pow((ax1 - ax2), 2))) < p) && ((sqrt(pow((ax1 - ax3), 2))) < p) && ((sqrt(pow((ax3 - ax2), 2))) < p)) && (((sqrt(pow((ay1 - ay2), 2))) < p) && ((sqrt(pow((ay1 - ay3), 2))) < p) && ((sqrt(pow((ay3 - ay2), 2))) < p))) return true;
    else return false;
}

bool dd_t1_t2(float M[6][2], float p) {
    float a1 = sqrt(pow((M[0][0] - M[3][0]), 2) + pow((M[0][1] - M[3][1]), 2));
    float a2 = sqrt(pow((M[2][0] - M[5][0]), 2) + pow((M[2][1] - M[5][1]), 2));
    float a3 = sqrt(pow((M[1][0] - M[4][0]), 2) + pow((M[1][1] - M[4][1]), 2));
    if ((a1 == a2 && a1 == a3) || (((sqrt(pow((a1 - a2), 2))) < p) && ((sqrt(pow((a1 - a3), 2))) < p) && ((sqrt(pow((a3 - a2), 2))) < p))) return true;
    else return false;
}

bool prr(float S[6][2], float p, ofstream& log) {
    if (dd_t1_t2(S, p)) {
        if (c_t1_t2(S, p)) {
            if (ds_t1_t2(S, p)) {
                if (v_t1_t2(S, p)) {
                    for (int i = 0; i < 6; i++) {
                        //log << S[i][0] << " " << S[i][1] << "\n";
                    }
                    log << "\n---\n";
                    return 1;
                }
                else {
                    //cout << "\n!4";
                    return 0;
                }
            }
            else {
                //cout << "\n!3";
                return 0;
            }
        }
        else {
            //cout << "\n!2";
            return 0;
        }
    }
    else {
        //cout << "\n!1";
        return 0;
    }
}

int poisk_6(float** mass, int n, ofstream& log, float p, float*** she, ofstream& out) {
    log << "\nЗаход в poisk_6\n";
    int q = 0;
    log << "Считано " << n << " точек\n";
    log << "Вот точки, которые могут быть в роли шестиугольника:\n";
    for (int i0 = 0; i0 < n; i0++) {
        for (int i1 = 0; i1 < n; i1++) {
            if (i0 == i1) continue;
            for (int i2 = 0; i2 < n; i2++) {
                if (i1 == i2) continue;
                for (int i3 = 0; i3 < n; i3++) {
                    if (i2 == i3) continue;
                    for (int i4 = 0; i4 < n; i4++) {
                        if (i3 == i4) continue;
                        for (int i5 = 0; i5 < n; i5++) {
                            if (i4 == i5) continue;
                            bool flag = false;

                            float SH[6][2] = { {mass[i0][0], mass[i0][1]}, {mass[i1][0], mass[i1][1]}, {mass[i2][0], mass[i2][1]}, {mass[i3][0], mass[i3][1]}, {mass[i4][0], mass[i4][1]}, {mass[i5][0], mass[i5][1]} };

                            float t1, t2, t3, t4;

                            // проверка на равенство точек
                            //for (int k = 0; k < 6; k++) {
                            //    for (int k2 = 0; k2 < 6; k2++) {
                            //        if (k == k2) break;
                            //        t1 = SH[k][0];
                            //        t2 = SH[k2][0];

                            //        if (t1 == t2) {
                            //            t3 = SH[k][1];
                            //            t4 = SH[k2][1];
                            //            if (t3 == t4) {
                            //                /*
                            //                log << "\nНекоторые точки совпадают!\n";
                            //                for (int f = 0; f < 6; f++) {
                            //                    log << "\n" << SH[f][0] << " " << SH[f][1];
                            //                }
                            //                log << "\n!!!!!!!!"; */
                            //                flag = true;
                            //            }
                            //        }
                            //    }
                            //}
                            //if (flag) break;

                            for (int f = 0; f < 6; f++) {
                                log << "\n" << SH[f][0] << " " << SH[f][1];
                            }

                            log << "\n----";

                            if (prr(SH, p, log)) {
                                log << "\nШестерка подошла! Записывание шестерки в массив правильных шестерок \
                                    -----------------------------\n";
                                for (int f = 0; f < 6; f++) {
                                    she[q][f][0] = SH[f][0];
                                    she[q][f][1] = SH[f][1];
                                }
                                q++;

                            }

                        }
                    }
                }
            }
        }
    }
    if (q > 0) {
        log << "\nПравильные шестерки: \n";
        out << "\nПравильные шестерки: \n";

        groupCyclicPermutations(she, q);
        for (unsigned i = 0; i < q / 6; i++) {
            for (unsigned j = 0; j < 6; j++) {
                out << she[i][j][0] << " " << she[i][j][1] <<  endl;
                log << she[i][j][0] << " " << she[i][j][1] << endl;
            }
            out << "..." <<  endl;
            log << "..." << endl;
        }


        //for (int i = 0; i < q; i++) {
        //    log << "Шестерка N " << i << ":" << endl;
        //    for (int f = 0; f < 6; f++) {
        //        log << she[i][f][0] << " " << she[i][f][1] << "\n";
        //        out << she[i][f][0] << " " << she[i][f][1] << "\n";
        //    }
        //    log << "...\n";
        //    out << "...\n";
        //}
    }
    return q;
}

int num_read(ifstream& in, ofstream& out, ofstream& log, float** mass, const int ogr) {
    log << "Заход в num_read\n";
    int count_n = 0, nn = 0;
    float n;
    char ch;
    log << "Координаты точек:\n";
    out << "Координаты точек:\n";
    while (count_n < ogr && !in.eof()) {
        nn = 0;
        while (nn < 2 && !in.eof()) {
            if (in.get(ch)) {
                if (ch == '\n') {
                    continue;
                }
                else {
                    in.seekg(-1, ios::cur);
                }
            }

            if (in >> n) {
                log << n << " ";
                out << n << " ";
                mass[count_n][nn] = n;
                nn++;
            }
            if (in.clear(), in.get() == '\n') continue;
        }
        log << "\n";
        out << "\n";
        count_n++;
    }
    if (count_n >= ogr) cout << "\nВ файле больше " << ogr << " координат :(\n";
    log << "-----------------------------------\n";
    return count_n;
}

int nc(ifstream& in) {
    int count_n = 0, nn = 0;
    float n;
    while (!in.eof()) {
        nn = 0;
        while (nn < 2 && !in.eof()) {
            if (in >> n) {
                nn++;
            }
            if (in.clear(), in.get() == '\n') continue;

        }
        count_n++;
    }
    return count_n;
}

int main() {
    //SetConsoleCP(1251);
    //SetConsoleOutputCP(1251);

    setlocale(LC_ALL, "rus");

    const string PATH1 = "C:\\Users\\Анечка\\Documents\\in.txt";
    const string PATH2 = "C:\\Users\\Анечка\\Documents\\out.txt";
    const string PATH3 = "C:\\Users\\Анечка\\Documents\\log.txt";

    ifstream points(PATH1, ios::in);
    ofstream out(PATH2, ios::out);
    ofstream log(PATH3, ios::out);

    if (!points) {
        cout << "Ошибка открытия входного файла.";
        return 0;
    }
    if (!out) {
        cout << "Ошибка открытия выходного файла.";
        return 0;
    }
    if (!log) {
        cout << "Ошибка открытия файла проверки.";
        return 0;
    }

    points.seekg(0, ios::end);
    if (points.tellg() == 0) {
        cout << "Ошибка (пустой файл).";
        out << "Ошибка (пустой файл).";
        return 0;
    }
    points.seekg(0, ios::beg);
    //-----------------------------------------------------------
    unsigned size = nc(points);
    float** ppp = new float* [size];
    for (int i = 0; i < size; i++) {
        ppp[i] = new float[2];
    }

    const int ogr = 50;

    points.close();
    points.open(PATH1);

    int coord_num = num_read(points, out, log, ppp, ogr);
    if (coord_num < 6) {
        log << "Количество точек меньше 6! Шестиугольник невозможен";
        return 0;
    }

    points.close();

    //!!!!!!!!!!!!!!!!!!!!!!
    float p = 0.0001;

    int rr = coord_num * coord_num;
    float*** she = new float** [rr] {};
    for (int i = 0; i < rr; i++) {
        she[i] = new float* [6] {};
        for (int ii = 0; ii < 6; ii++) {
            she[i][ii] = new float[2] {};
        }
    }

    //q - кол-во шестерок
    int q = poisk_6(ppp, coord_num, log, p, she, out) / 6;

    float*** sorted_she = new float** [rr / 6];
    for (unsigned i = 0; i < rr / 6; i++)
    {
        sorted_she[i] = new float* [6];
        for (unsigned j = 0; j < 6; j++)
        {
            sorted_she[i][j] = new float[2];
            sorted_she[i][j][0] = she[i * 6][j][0];
            sorted_she[i][j][1] = she[i * 6][j][1];
        }
    }

    if (q == 0) {
        log << "\nПравильные шестиугольники не найдены";
        out << "\nПравильные шестиугольники не найдены";
        return 0;
    }

    if (q == 1 && coord_num == 6) {
        log << "\nВ файле содержится один шестиугольник, не содержащий точек";
        out << "\nВ файле содержится один шестиугольник, не содержащий точек";
        return 0;
    }

    float* p_pro = new float [q] {};

    log.close();
    ofstream logh(PATH3, ios::app);

    kal(ppp, p_pro, sorted_she, q, coord_num, p, logh);

    logh.close();
    ofstream logq(PATH3, ios::app);

    float mx = mmm(p_pro, q);
    pr_m(sorted_she, logq, mx, out, p_pro, q);

    //--------------------------------------------
    logq.close();
    out.close();

    for (int i = 0; i < ogr; i++) {
        delete[] ppp[i];
    }
    delete[] ppp;
    delete[] p_pro;

    /*for (int i = 0; i < q; i++) {
        for (int ii = 0; ii < 6; ii++) {
            delete[] she[i][ii];
        }
        delete[] she[i];
    }*/

    return 0;
}