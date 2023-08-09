#include <string>
#include <fstream>
#include <vector>

using namespace std;
class FileHandler
{
public:


    vector<string> readbpRNAFile(string filename)
    {
        string line;
        string bpRNAseq, DB;
        ifstream myfile(filename);
        if (myfile.is_open())
        {
            int counter = 0;
            while (getline(myfile, line))
            {
                if (counter == 0)
                {
                    bpRNAseq = line;
                }
                if (counter == 1)
                {
                    DB = line;
                }
                counter++;
            }
        }

        return { bpRNAseq, DB };
    }
};