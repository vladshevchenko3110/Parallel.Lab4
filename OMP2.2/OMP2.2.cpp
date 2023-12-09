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
    using second_t = chrono::duration<double, ratio<1, 1000> >;
    chrono::time_point<clock_t> m_beg;
public:
    Timer() : m_beg(clock_t::now()) {}
    void reset() { m_beg = clock_t::now(); }
    double elapsed() const {
        return chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
    }
};

void SearchComparison(int loadarraysize) {

    cout << "======== Comparison with an array of " << loadarraysize << " elements ========" << endl;
    vector<double> loadedArray(loadarraysize);
    ifstream inFile("array.bin", ios::binary);
    inFile.read(reinterpret_cast<char*>(loadedArray.data()), loadarraysize * sizeof(double));
    inFile.close();

    
    Timer t1;
    int min = loadedArray[0]; 
    for (int i = 1; i < loadarraysize; ++i) {
        if (loadedArray[i] < min) {
            min = loadedArray[i]; 
        }
    }
    double time1 = t1.elapsed();
    cout << "Sequential processing time : " << time1 << " ms" << '\n';
    cout << "Minimum array value: " << min << endl;
    

    ofstream newOutFile("serial_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream newOutFileMinResult("serial_new_minresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    
    newOutFile.write(reinterpret_cast<char*>(loadedArray.data()), loadarraysize * sizeof(double));

    newOutFile.close();

    newOutFileMinResult.write(reinterpret_cast<char*>(&min), sizeof(double));
    newOutFileMinResult.close();
    


    vector<double> omploadedArray(loadarraysize);
    ifstream OMPinFile("array.bin", ios::binary);
    OMPinFile.read(reinterpret_cast<char*>(omploadedArray.data()), loadarraysize * sizeof(double));
    OMPinFile.close();

   
    int ompminvalue = omploadedArray[0];
    int localmin = omploadedArray[0]; 
    
    
    Timer t2;
    #pragma omp parallel 
    {     

        #pragma omp for 
        for (int i = 0; i < loadarraysize; ++i) {
            if (omploadedArray[i] < localmin) {
                localmin = omploadedArray[i];
            }
        }
        
        #pragma omp critical
        {
            if (localmin < ompminvalue) {
                ompminvalue = localmin;
            }
        }
    }
    double time2 = t2.elapsed();

    cout << "Parallel processing time : " << time2 << " ms" << '\n';
    cout << "Minimum array value: " << ompminvalue << endl;
    
    cout << "Acceleration of parallel block processing: " << time1 / time2 << endl;
    ofstream OMPnewOutFile("parallel_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFileMinResult("parallel_new_minresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    

    OMPnewOutFile.write(reinterpret_cast<char*>(omploadedArray.data()), loadarraysize * sizeof(double));

    OMPnewOutFile.close();

    OMPnewOutFileMinResult.write(reinterpret_cast<char*>(&ompminvalue), sizeof(double));
    OMPnewOutFileMinResult.close();

    vector<double> omploadedArray1(loadarraysize);
    ifstream OMPinFile1("array.bin", ios::binary);
    OMPinFile1.read(reinterpret_cast<char*>(omploadedArray1.data()), loadarraysize * sizeof(double));
    OMPinFile1.close();


    int ompminvalue1 = omploadedArray1[0];
    int localmin1 = omploadedArray1[0];
    
    Timer t3;
    #pragma omp parallel for reduction(min : ompminvalue1)
    for (int i = 0; i < loadarraysize; ++i) {        
        if (omploadedArray1[i] < ompminvalue1) {
            ompminvalue1 = omploadedArray1[i];
        }
    }
    double time3 = t3.elapsed();

    cout << "Parallel processing time 2: " << time3 << " ms" << '\n';
    cout << "Minimum array value: " << ompminvalue1 << endl;

    cout << "Acceleration of parallel block processing 2: " << time1 / time3 << endl;
    ofstream OMPnewOutFile1("parallel_new_array1_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFileMinResult1("parallel_new_minresult1_" + to_string(loadarraysize) + ".bin", ios::binary);


    OMPnewOutFile1.write(reinterpret_cast<char*>(omploadedArray.data()), loadarraysize * sizeof(double));
    
    OMPnewOutFile1.close();

    OMPnewOutFileMinResult1.write(reinterpret_cast<char*>(&ompminvalue), sizeof(double));
    OMPnewOutFileMinResult1.close();
   
    cout << endl;
}


int main(int argc, char* argv[])
{
    /*srand((unsigned int)time(0));


    const int arraySize = 100000000;


    vector<double> array(arraySize);
    for (int i = 0; i < arraySize; ++i) {
        array[i] = ((rand() % 2000001)-1000000);
    }


    ofstream outFile("array.bin", ios::binary);
    outFile.write(reinterpret_cast<char*>(array.data()), arraySize * sizeof(double));
    outFile.close();*/

    SearchComparison(stoi(argv[1]));
    SearchComparison(stoi(argv[2]));
    SearchComparison(stoi(argv[3]));
    SearchComparison(stoi(argv[4]));
    SearchComparison(stoi(argv[5]));
    SearchComparison(stoi(argv[6]));
    SearchComparison(stoi(argv[7]));

}