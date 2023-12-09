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

void CountingComparison(int loadarraysize) {

    cout << "======== Comparison with an array of " << loadarraysize << " elements ========" << endl;
    vector<double> loadedArray(loadarraysize);
    ifstream inFile("array.bin", ios::binary);    
    inFile.read(reinterpret_cast<char*>(loadedArray.data()), loadarraysize * sizeof(double));
    inFile.close();

    int poscount = 0;
    int negcount = 0;
    int nullcount = 0;
    Timer t1;
    for (int i = 0; i < loadarraysize; ++i) {      
        
        if (loadedArray[i] > 0) {            
            poscount++;
        } 
        else if (loadedArray[i] < 0) {
            negcount++;
        } 
        else {            
            nullcount++;
        }    
    }
    double time1 = t1.elapsed();
    cout << "Sequential processing time : " << time1 << " ms" << '\n';
    cout << "Number of positive values: " << poscount << endl;
    cout << "Number of negative values: " << negcount << endl;
    cout << "Number of zero values: " << nullcount << endl;

    ofstream newOutFile("serial_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream newOutFilePosResult("serial_new_posresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream newOutFileNegResult("serial_new_negresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream newOutFileNullResult("serial_new_nullresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    
    newOutFile.write(reinterpret_cast<char*>(loadedArray.data()), loadarraysize * sizeof(double));
    
    newOutFile.close();

    newOutFilePosResult.write(reinterpret_cast<char*>(&poscount), sizeof(int));
    newOutFilePosResult.close();
    newOutFileNegResult.write(reinterpret_cast<char*>(&negcount), sizeof(int));
    newOutFileNegResult.close();
    newOutFileNullResult.write(reinterpret_cast<char*>(&nullcount), sizeof(int));
    newOutFileNullResult.close();


    vector<double> omploadedArray(loadarraysize);
    ifstream OMPinFile("array.bin", ios::binary);
    OMPinFile.read(reinterpret_cast<char*>(omploadedArray.data()), loadarraysize * sizeof(double));
    OMPinFile.close();

    int ompposcount = 0;
    int ompnegcount = 0;
    int ompnullcount = 0;
   
    
    Timer t2;      
    #pragma omp parallel for 
    for (int i = 0; i < loadarraysize; i++) {
        
        if (omploadedArray[i] > 0) {
            #pragma omp atomic
            ompposcount++;  
        }       
        else if (omploadedArray[i] < 0) {
            #pragma omp atomic
            ompnegcount++;
        }           
        else {
            #pragma omp atomic
            ompnullcount++;  
        }              
               
    }        
    double time2 = t2.elapsed();
    cout << "Parallel processing time : " << time2 << " ms" << '\n';
    cout << "Number of positive values: " << ompposcount << endl;
    cout << "Number of negative values: " << ompnegcount << endl;
    cout << "Number of zero values: " << ompnullcount << endl;
    cout << "Acceleration of parallel block processing: " << time1 / time2 << endl;
    ofstream OMPnewOutFile("parallel_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFilePosResult("parallel_new_posresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFileNegResult("parallel_new_negresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFileNullResult("parallel_new_nullresult_" + to_string(loadarraysize) + ".bin", ios::binary);

    OMPnewOutFile.write(reinterpret_cast<char*>(omploadedArray.data()), loadarraysize * sizeof(double));

    OMPnewOutFile.close();

    OMPnewOutFilePosResult.write(reinterpret_cast<char*>(&ompposcount), sizeof(int));
    OMPnewOutFilePosResult.close();
    OMPnewOutFileNegResult.write(reinterpret_cast<char*>(&ompnegcount), sizeof(int));
    OMPnewOutFileNegResult.close();
    OMPnewOutFileNullResult.write(reinterpret_cast<char*>(&ompnullcount), sizeof(int));
    OMPnewOutFileNullResult.close();

    vector<double> omploadedArray1(loadarraysize);
    ifstream OMPinFile1("array.bin", ios::binary);
    OMPinFile1.read(reinterpret_cast<char*>(omploadedArray1.data()), loadarraysize * sizeof(double));
    OMPinFile1.close();

    int ompposcount1 = 0;
    int ompnegcount1 = 0;
    int ompnullcount1 = 0;

    
    Timer t3;
    #pragma omp parallel for reduction(+:ompposcount1,ompnegcount1,ompnullcount1)
    for (int i = 0; i < loadarraysize; i++) {

        if (omploadedArray1[i] > 0) {
            ompposcount1++;
        }
        else if (omploadedArray1[i] < 0) {
            ompnegcount1++;
        }
        else {
            ompnullcount1++;
        }

    }
    double time3 = t3.elapsed();
    cout << "Parallel processing time : " << time3 << " ms" << '\n';
    cout << "Number of positive values: " << ompposcount1 << endl;
    cout << "Number of negative values: " << ompnegcount1 << endl;
    cout << "Number of zero values: " << ompnullcount1 << endl;
    cout << "Acceleration of parallel block processing 2: " << time1 / time3 << endl;
    ofstream OMPnewOutFile1("parallel_new_array_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFilePosResult1("parallel_new_posresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFileNegResult1("parallel_new_negresult_" + to_string(loadarraysize) + ".bin", ios::binary);
    ofstream OMPnewOutFileNullResult1("parallel_new_nullresult_" + to_string(loadarraysize) + ".bin", ios::binary);

    OMPnewOutFile1.write(reinterpret_cast<char*>(omploadedArray1.data()), loadarraysize * sizeof(double));

    OMPnewOutFile1.close();

    OMPnewOutFilePosResult1.write(reinterpret_cast<char*>(&ompposcount1), sizeof(int));
    OMPnewOutFilePosResult1.close();
    OMPnewOutFileNegResult1.write(reinterpret_cast<char*>(&ompnegcount1), sizeof(int));
    OMPnewOutFileNegResult1.close();
    OMPnewOutFileNullResult1.write(reinterpret_cast<char*>(&ompnullcount1), sizeof(int));
    OMPnewOutFileNullResult1.close();

    cout << endl;
}


int main(int argc, char* argv[])
{
    /*srand((unsigned int)time(0));


    const int arraySize = 100000000;


    vector<double> array(arraySize);
    for (int i = 0; i < arraySize; ++i) {
        array[i] = ((rand() % 20001)-10000);
    }


    ofstream outFile("array.bin", ios::binary);
    outFile.write(reinterpret_cast<char*>(array.data()), arraySize * sizeof(double));    
    outFile.close();*/

    
    CountingComparison(stoi(argv[1]));
    CountingComparison(stoi(argv[2]));
    CountingComparison(stoi(argv[3]));
    CountingComparison(stoi(argv[4]));
    CountingComparison(stoi(argv[5]));
    CountingComparison(stoi(argv[6]));
    CountingComparison(stoi(argv[7]));

}