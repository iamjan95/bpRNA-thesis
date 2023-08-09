#include <string>
#include <fstream>
#include <vector>
#include <map>
using namespace std;

class FileInputs
{
public:
    void writetoFile(string bpRNAseq, string filename,string DB)
    {
        ofstream myfile(filename);
        if (myfile.is_open())
        {
            myfile << bpRNAseq<<endl<<DB;
        }
    }

    void readdbnFile(string filename)
    {
        string line;
        ifstream myfile(filename);
        if (myfile.is_open())
        {
            int counter = 0;
            while (getline(myfile, line))
            {
                if (counter == 4)
                {
                    DB = line;
                }
                counter++;
            }
        }
        myfile.close();
        fixDB();
    }

    void fixDB()
    {

        for (int i = 0; i < DB.length(); i++)
        {
            if (DB[i] != ')' && DB[i] != '(')
                DB[i] = '.';
        }
    }

    void readBPSeqFile(string filename)
    {
        string line;
        ifstream myfile(filename);
        string delimiter = " ";
        if (myfile.is_open())
        {

            while (getline(myfile, line))
            {
                string token;
                int seqnum;
                int pos = 0;
                int counter = 0;
                while ((pos = line.find(delimiter)) != std::string::npos)
                {
                    token = line.substr(0, pos);
                    if (counter == 0)
                        seqnum = stoi(token);
                    if (counter == 1)
                        seq.push_back(token[0]);
                    if (counter == 2)
                    {
                        bp[seqnum] = stoi(token);
                        bp[stoi(token)] = seqnum;
                    }
                    line.erase(0, pos + delimiter.length());
                    counter++;
                }
            }
            myfile.close();
        }
    }

    string getDBinfo()
    {
        return DB;
    }
    
    vector<char> getseq()
    {
        return seq;
    }

    map <int, int> getBPinfo()
    {
        return bp;
    }

  

private:
    string DB;
    vector<char> seq;
    map <int, int> bp;
};