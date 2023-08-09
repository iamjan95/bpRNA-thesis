#include "bpRNAsequencer.cpp"
#include<iostream>

int main(int argc, char** argv)
{
    FileInputs fileHandler;
    bpRNASequencer sequencer;
    if (argc > 2)
    {
        
        string filename = argv[1];
        if (filename.substr(filename.find_last_of(".") + 1) == "dbn") 
        {
            fileHandler.readdbnFile(filename);
            sequencer.setDB(fileHandler.getDBinfo());
            sequencer.convertDB();
        }
        if (filename.substr(filename.find_last_of(".") + 1) == "bpseq")
        {
            fileHandler.readBPSeqFile(filename);
            sequencer.setDB(fileHandler.getDBinfo());
            sequencer.setseq(fileHandler.getseq());
            sequencer.convertbpSeq();
        }
        string bpRNAfilename = argv[2];
        fileHandler.writetoFile(sequencer.getStructureVector(), bpRNAfilename,sequencer.getDB());
    }
}