
#include "Align.cpp"
#include "FileHandler.cpp"

string editbpRNAseq(string bpRNAseq, string dotbracket)
{
    string seq = "";
    string right_symbols = ")>]}";
    string left_symbols = "([<{";
    for (int i = 0; i < bpRNAseq.length(); i++) {
        if ((right_symbols.find(dotbracket[i]) != string::npos || islower(dotbracket[i])) && bpRNAseq[i] == 'S') {
            seq += 'R';
        }
        else if ((left_symbols.find(dotbracket[i]) != string::npos || isupper(dotbracket[i])) && bpRNAseq[i] == 'S') {
            seq += 'L';
        }
        else
        {
            seq += bpRNAseq[i];
        }
    }
    return seq;
}

int main(int argc, char** argv)
{
        FileHandler fileH;
        int bandwidth;
        if (argc > 2)
        {

            string filename1 = argv[1];
            string filename2 = argv[2];
            if (argc == 4)
                bandwidth = stoi(argv[3]);
            else
                bandwidth = 20;
            vector<string> file1Info = fileH.readbpRNAFile(filename1);
            vector<string> file2Info = fileH.readbpRNAFile(filename2);
            string bpRNA1 = file1Info[0];
            string db1 = file1Info[1];
            bpRNA1 = editbpRNAseq(bpRNA1, db1);
            string bpRNA2 = file2Info[0];
            string db2 = file2Info[1];
            bpRNA2 = editbpRNAseq(bpRNA2, db2);
            int row = bpRNA1.length() + 1;
            int col = bpRNA2.length() + 1;

            Align a(row, col, bpRNA1, bpRNA2);
            a.setbandwidth(bandwidth);
            a.performAlignment();
            a.traceback();
            cout << endl << "Distance: " << a.calcDistance();
        }
}



