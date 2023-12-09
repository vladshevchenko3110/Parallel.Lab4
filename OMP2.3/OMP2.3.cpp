#include <iostream>
#include <omp.h>
#include <time.h>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

class Timer {
private:
    using clock_t = chrono::high_resolution_clock;
    using second_t = chrono::duration<double, ratio<1,1000> >;
    chrono::time_point<clock_t> m_beg;
public:
    Timer() : m_beg(clock_t::now()) {}
    void reset() { m_beg = clock_t::now(); }
    double elapsed() const {
        return chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }
};

void bubbleSort(vector<double>& array) {
    int size = array.size();

    int phase, i, temp;
    for (phase = 0; phase < size; ++phase)
    {
        if (phase % 2 == 0) 
        {
            for (i = 1; i < size; i += 2)
                if (array[i - 1] > array[i])
                {
                    temp = array[i];
                    array[i] = array[i - 1];
                    array[i - 1] = temp;
                }
        }
        else 
        {
            for (i = 1; i < size - 1; i += 2)
                if (array[i] > array[i + 1])
                {
                    temp = array[i];
                    array[i] = array[i + 1];
                    array[i + 1] = temp;
                }
        }
    }

}


void parallelBubbleSort(vector<double>& array) {

    int size = array.size();
    int phase, i, temp;
    for (phase = 0; phase < size; ++phase)
    {
        if (phase % 2 == 0) 
        {
            #pragma omp parallel for shared(array,size) private(i,temp)
            for (i = 1; i < size; i += 2)
                if (array[i - 1] > array[i])
                {
                    temp = array[i];
                    array[i] = array[i - 1];
                    array[i - 1] = temp;
                }
        }
        else 
        {
            #pragma omp parallel for shared(array,size) private(i,temp)
            for (i = 1; i < size - 1; i += 2)
                if (array[i] > array[i + 1])
                {
                    temp = array[i];
                    array[i] = array[i + 1];
                    array[i + 1] = temp;
                }
        }
    }

}

int partition(vector<double>& array, int low, int high) {
    double pivot = array[high];
    int i = (low - 1);

    for (int j = low; j <= high - 1; j++) {
        if (array[j] < pivot) {
            i++;
            swap(array[i], array[j]);
        }
    }

    swap(array[i + 1], array[high]);
    return (i + 1);
}

void quickSort(vector<double>& array, int low, int high) {
    if (low < high) {
        int pi = partition(array, low, high);

        quickSort(array, low, pi - 1);
        quickSort(array, pi + 1, high);
    }
}

void parallelQuickSort(vector<double>& array, int low, int high, int low_limit, int threads) {
    if (low < high) {
        
        if ((high - low) < low_limit || threads == 1) {
            quickSort(array, low, high);
        }
        else {
            int pi = partition(array, low, high);
            #pragma omp parallel sections 
            {
                
                #pragma omp section
                parallelQuickSort(array, low, pi - 1, low_limit, threads / 2);
                #pragma omp section
                parallelQuickSort(array, pi + 1, high, low_limit, threads - threads / 2);
            }
        }        
                        
    }
}

void SortComparison(int loadarraysize) {

    
    cout << "======== Comparison with an array of " << loadarraysize << " elements ========" << endl;

    vector<double> loadedArray1(loadarraysize);
    ifstream inFile1("array.bin", ios::binary);
    inFile1.read(reinterpret_cast<char*>(loadedArray1.data()), loadarraysize * sizeof(double));
    inFile1.close();

    vector<double> loadedArray2(loadarraysize);
    ifstream inFile2("array.bin", ios::binary);
    inFile2.read(reinterpret_cast<char*>(loadedArray2.data()), loadarraysize * sizeof(double));
    inFile2.close();

    double time1 = 0;
    Timer t1;
    if (loadarraysize < 100000) {
        bubbleSort(loadedArray1);
        time1 = t1.elapsed();
        cout << "Sequential bubble sort time : " << time1 << " ms" << '\n';
        
    }
    t1.reset();
    quickSort(loadedArray2, 0, loadarraysize - 1);
    double time3 = t1.elapsed();
    cout << "Sequential quick sort time : " << time3 << " ms" << '\n';


    ofstream newOutFile1("serial_bubble_sort_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream newOutFile2("serial_quick_sort_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);

    newOutFile1.write(reinterpret_cast<char*>(loadedArray1.data()), loadarraysize * sizeof(double));
    newOutFile1.close();

    newOutFile2.write(reinterpret_cast<char*>(loadedArray1.data()), loadarraysize * sizeof(double));
    newOutFile2.close();    



    vector<double> omploadedArray1(loadarraysize);
    ifstream OMPinFile1("array.bin", ios::binary);
    OMPinFile1.read(reinterpret_cast<char*>(omploadedArray1.data()), loadarraysize * sizeof(double));
    OMPinFile1.close();

    vector<double> omploadedArray2(loadarraysize);
    ifstream OMPinFile2("array.bin", ios::binary);
    OMPinFile2.read(reinterpret_cast<char*>(omploadedArray2.data()), loadarraysize * sizeof(double));
    OMPinFile2.close();


    double time2 = 0;
    
    Timer t2;
    if (loadarraysize < 100000) {
        
        parallelBubbleSort(omploadedArray1);
                
        time2 = t2.elapsed();
        cout << "Parallel bubble sort time : " << time2 << " ms" << '\n';        
    }    
    t2.reset();     
    parallelQuickSort(omploadedArray2, 0, loadarraysize - 1,1000,12);
      
    double time4 = t2.elapsed();
    cout << "Parallel quick sort time : " << time4 << " ms" << '\n';

    if (loadarraysize < 100000) {
        cout << "Acceleration of parallel block bubble sort: " << time1 / time2 << endl;
    }    
    cout << "Acceleration of parallel block quick sort: " << time3 / time4 << endl;

    ofstream OMPnewOutFile1("parallel_bubble_sort_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFile2("parallel_quick_sort_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);


    OMPnewOutFile1.write(reinterpret_cast<char*>(omploadedArray1.data()), loadarraysize * sizeof(double));
    OMPnewOutFile1.close();
    OMPnewOutFile2.write(reinterpret_cast<char*>(omploadedArray2.data()), loadarraysize * sizeof(double));
    OMPnewOutFile2.close();

    

    cout << endl;
}


int main(int argc, char* argv[])
{
    
    omp_set_nested(1);
    /*srand((unsigned int)time(0));
    

    const int arraySize = 100000000;


    vector<double> array(arraySize);
    for (int i = 0; i < arraySize; ++i) {
        array[i] = rand() % 100000;
    }


    ofstream outFile("array.bin", ios::binary);
    outFile.write(reinterpret_cast<char*>(array.data()), arraySize * sizeof(double));
    outFile.close();*/

    SortComparison(stoi(argv[1]));
    SortComparison(stoi(argv[2]));
    SortComparison(stoi(argv[3]));
    SortComparison(stoi(argv[4]));
    SortComparison(stoi(argv[5]));
    SortComparison(stoi(argv[6]));
    

}