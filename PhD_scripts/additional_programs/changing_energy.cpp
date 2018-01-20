	

    #include <iostream>
    #include <fstream>
    #include <cstdlib>
    #include <vector>
    #include <iomanip>
    using namespace std;
     
     
    double extractFromBuffor(char* buffor, int column){
     
      column--;
     
        while(column > 0){
            if(buffor[0] == ' '){
                column--;
                while(buffor[0] == ' '){
                    buffor++;
                }
            }
            buffor++;
        }
     
        return atof(buffor);
     
    }
     
    void replaceInBuffor(char* buffor, int column, double data){
     
        column--;
     
        while(column > 0){
            if(*buffor == ' '){
                column--;
            }
            buffor++;
        }
     
     
     
    }
     
     
    int main()
    {
     
        ifstream inFile;
        string inFileName;
        inFileName ="prop.out";
        ofstream outFile;
     
     
        char buffor[1024];
     
        vector<double> data1;
        vector<double> data5;
     
        double foo = 0;
     
        inFile.open(inFileName.c_str());
        outFile.open("changing.prop.out");
     
        int i1 = 0;
        int i5 = 0;
        while(!inFile.eof()){
            inFile.getline(buffor,1024);
     
        for(int j = 0; j<1024; j++){
           if (buffor[j] == 'D'){
               buffor[j] = 'e';
           }
        }

            switch(buffor[0]){
                case '5':{
                    data5.push_back(extractFromBuffor(buffor,9));
                    foo = data5[i5] - data5[0];
                    foo *=  27.211;
                    i5++;
                    break;
                }
                case '1':{
                    data1.push_back(extractFromBuffor(buffor,9));
                    foo = data1[i1];
                    foo /=  0.3934;
                    i1++;
                    break;
                }
     
                default:
                    break;
     
            }
     
     
        outFile << buffor<<"  -  " <<scientific <<foo<< endl;
     
        }
     
     
     
        //system("pause");
    }


