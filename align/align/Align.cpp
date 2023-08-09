#include "Constants.h"

class Align
{
public:
    
    double generate_score(int i, int j)
    {
        char a = seq1[i - 1];
        char b = seq2[j - 1];
        if (score.find(make_pair(a, b)) != score.end())
            return score[(make_pair(a, b))];
        else
            return score[(make_pair(b, a))];
    }

    bool inside_band(int i, int j)
    {
        int diag = int(ceil(((float(rows) - 1) / (float(cols) - 1)) * j));
        int rangemin = diag - bandwidth;
        int rangemax = diag + bandwidth;
        if (i >= rangemin && i <= rangemax)
            return true;
        else
            return false;
    }

    char getDirectionChar(int x)
    {
        switch (x)
        {
        case 0: return 'm';
        case 1: return 'x';
        case 2: return 'y';
        }
    }
    int getDirectionInt(char x)
    {
        switch (x)
        {
        case 'm': return 0;
        case 'x': return 1;
        case 'y': return 2;

        }
    }


    void performAlignment()
    {
        int count = -1;
        for (int i = 1; i <= bandwidth; i++)
        {
            middleDir[i][0] = 'y';
            xDir[i][0] = 'y';
            yDir[i][0] = 'y';
            xMat[i][0] = GAP + (i - 1) * EXTEND;

        }
        for (int j = 1; j <= bandwidth; j++)
        {
            middleDir[0][j] = 'x';
            xDir[0][j] = 'x';
            yDir[0][j] = 'x';
            yMat[0][j] = GAP + (j - 1) * EXTEND;
        }
        for (int j = 1; j < cols; j++)
        {
            for (int d = -bandwidth; d < bandwidth; d++)
            {
                int i = int(ceil(((float(rows) - 1) / (cols - 1)) * j + d));
                if (i >= 1 && i < rows)
                {
                    calulateScore(i, j);
                    middleMat[i][j] = score_m[0];
                    xMat[i][j] = score_x[0];
                    yMat[i][j] = score_y[0];
                    middleDir[i][j] = getDirectionChar(score_m[1]);
                    xDir[i][j] = getDirectionChar(score_x[1]);
                    yDir[i][j] = getDirectionChar(score_y[1]);
                }
            }
        }
        findFinalmax();
    }

    void setbandwidth(int bw)
    {
        bandwidth = bw;
    }

    vector<double> getMaxPosition(vector<double> mxy) 
    {
        double max_score = -1e8;
        vector<int> max_index;
        bool flag_m = 0;

        for (int i = 0; i < 3; i++) {
            if (mxy[i] > max_score) {
                max_score = mxy[i];
            }
        }
        for (int i = 0; i < 3; i++)
        {
            if (abs(mxy[i] - max_score) < 0.0001)
            {
                if (i == 0)
                {
                    flag_m = 1;
                }
                max_index.push_back(i);
            }
        }
        if (flag_m == 1)
        {
            return { max_score, 0 };  
        }
        else if (max_index.size() > 1)
        {
            return { max_score, 2 }; 
        }
        else if (max_index.size() == 1 && max_index[0] == 1)
        {
            return { max_score, 1 }; 
        }
        else if (max_index.size() == 1 && max_index[0] == 2)
        {
            return { max_score, 2 };
        }
    }

    double calcDistance()
    {
        int coun = 0;
        for (int i = 0; i < Aseq1.length(); i++)
        {
            if (Aseq1[i] == Aseq2[i])
                coun++;
        }
        double dist = (double)coun / Aseq1.length();
        return 1 - dist;
    }


    void findFinalmax()
    {
        int middleMax = middleMat[rows - 1][cols - 1];
        int X_max = xMat[rows - 1][cols - 1];
        int Y_max = yMat[rows - 1][cols - 1];
        vector<vector<int>> list = { {X_max, 'x'}, {Y_max, 'y'}, {middleMax, 'm'} };
        int score = INT_MIN;
        maxDir = 'n';
        for (auto i : list)
        {
            int value = i[0];
            char matrix = i[1];
            if (value > score) {
                score = value;
                maxDir = matrix;
            }
            else if (value == score && maxDir != 'm') {
                if (matrix == 'm') {
                    score = value;
                    maxDir = matrix;
                }
            }
        }
    }
    void traceback()
    {
        Aseq1 = "";
        Aseq2 = "";
        //cout << middleMat[rows - 1][cols - 1] << endl;
        int i = rows - 1, j = cols - 1;
        char direction;
        char current_matrix;
        if (maxDir == 'm') 
        {
            direction = middleDir[i][j];
        }
        else if (maxDir == 'y') 
        {
            direction = yDir[i][j];
        }
        else if (maxDir == 'x') 
        {
            direction = xDir[i][j];
        }
        current_matrix = maxDir;
        while (!(i == 0 && j == 0)) {
            if (i == 0 && j > 0) {
                Aseq1 = '-' + Aseq1;
                Aseq2 = seq2[j - 1] + Aseq2;

            }
            else if (j == 0 && i > 0) {
                Aseq1 = seq1[i - 1] + Aseq1;
                Aseq2 = '-' + Aseq2;

            }
            else if (current_matrix == 'm') {
                Aseq1 = seq1[i - 1] + Aseq1;
                Aseq2 = seq1[j - 1] + Aseq2;

            }
            else if (current_matrix == 'y') {
                Aseq1 = seq1[i - 1] + Aseq1;
                Aseq2 = '-' + Aseq2;

            }
            else if (current_matrix == 'x') {
                Aseq1 = '-' + Aseq1;
                Aseq2 = seq2[j - 1] + Aseq2;

            }
            vector<int> nMove = nextmove(i, j, current_matrix, direction);
            i = nMove[0];
            j = nMove[1];
            current_matrix = getDirectionChar(nMove[2]);
            direction = getDirectionChar(nMove[3]);

        }


        cout << Aseq1 << endl;
        cout << Aseq2 << endl;
    }

    vector<int> nextmove(int i, int j, char current_matrix, char direction)
    {
        vector<int> nMove = { 0,0,0,0 };
        if (j == 0 && i > 0)
        {
            nMove = { i - 1,j,2,2 };
        }
        else if (i == 0 && j > 0)
        {
            nMove = { i,j - 1,1,1 };
        }
        else if (i > 0 && j > 0)
        {
            if (current_matrix == 'm')
            {
                nMove[2] = getDirectionInt(direction);
                nMove[0] = i - 1;
                nMove[1] = j - 1;
                if (direction == 'm')
                    nMove[3] = getDirectionInt(middleDir[i][j]);
                else if (direction == 'y')
                    nMove[3] = getDirectionInt(yDir[i][j]);
                else if (direction == 'x')
                    nMove[3] = getDirectionInt(xDir[i][j]);
            }
            else if (current_matrix == 'x')
            {
                nMove[2] = getDirectionInt(direction);
                nMove[0] = i;
                nMove[1] = j - 1;
                if (direction == 'm')
                    nMove[3] = getDirectionInt(middleDir[i][j]);
                else if (direction == 'y')
                    nMove[3] = getDirectionInt(yDir[i][j]);
                else if (direction == 'x')
                    nMove[3] = getDirectionInt(xDir[i][j]);
            }
            else if (current_matrix == 'y')
            {
                nMove[2] = getDirectionInt(direction);
                nMove[0] = i - 1;
                nMove[1] = j;
                if (direction == 'm')
                    nMove[3] = getDirectionInt(middleDir[i][j]);
                else if (direction == 'y')
                    nMove[3] = getDirectionInt(yDir[i][j]);
                else if (direction == 'x')
                    nMove[3] = getDirectionInt(xDir[i][j]);
            }
        }
        return nMove;
    }

    void calulateScore(int i, int j)
    {
        vector<double> mxy(3);

        double score = generate_score(i, j);
        if (inside_band(i - 1, j))
        {
            int extendup = getGapExtendScore(i, j, Direction::UP);
            int gapup = getGapScore(i, j, Direction::UP);
            mxy[2] = yMat[i - 1][j] + extendup;
            mxy[1] = xMat[i - 1][j] + gapup;
            mxy[0] = middleMat[i - 1][j] + gapup;

            score_y = getMaxPosition(mxy);
        }
        else if (!inside_band(i - 1, j))
        {
            double extendup = getGapExtendScore(i, j, Direction::UP);
            double score_yy = yMat[i - 1][j - 1] + extendup;
            score_y = { score_yy,2 };
        }
        if (inside_band(i, j - 1))
        {
            double extendleft = getGapExtendScore(i, j, Direction::LEFT);
            double gapleft = getGapScore(i, j, Direction::LEFT);
            mxy[2] = yMat[i][j - 1] + extendleft;
            mxy[1] = xMat[i][j - 1] + gapleft;
            mxy[0] = middleMat[i][j - 1] + gapleft;
            score_x = getMaxPosition(mxy);

        }
        else if (!inside_band(i, j - 1))
        {
            double extendleft = getGapExtendScore(i, j, Direction::LEFT);
            double score_xx = xMat[i - 1][j - 1] + extendleft;
            score_x = { score_xx,1 };
        }


        mxy[0] = middleMat[i - 1][j - 1] + score;
        mxy[1] = xMat[i - 1][j - 1] + score;
        mxy[2] = yMat[i - 1][j - 1] + score;
        score_m = getMaxPosition(mxy);
    }

    double getGapScore(int i, int j, Direction dir)
    {
        if (dir == Direction::UP)
        {
            char b = seq1[i - 1];
            char a = '-';
            if (gap_score.find(make_pair(a, b)) != gap_score.end() && i - 2 >= 0 && j - 2 >= 0 && i < seq1.length() && j < seq2.length())
            {
                if (seq2[j] == seq1[i] && seq1[i - 2] == seq2[j - 1] && seq1[i] == seq2[j - 1])
                    return gap_score[make_pair(a, b)];
                else
                    return GAP;
            }
            else
                return GAP;
        }
        else if (dir == Direction::LEFT)
        {
            char b = seq1[i - 1];
            char a = '-';
            if (gap_score.find(make_pair(a, b)) != gap_score.end() && i - 2 >= 0 && j - 2 >= 0 && i < seq1.length() && j < seq2.length())
            {
                if (seq2[j] == seq1[i] && seq1[i - 1] == seq2[j - 2] && seq1[i] == seq2[j - 1])
                    return gap_score[make_pair(a, b)];
                else
                    return GAP;
            }
            else
                return GAP;
        }

    }

    double getGapExtendScore(int i, int j, Direction dir)
    {
        if (dir == Direction::UP)
        {
            char b = seq1[i - 1];
            char a = '-';
            if (extend_score.find(make_pair(a, b)) != extend_score.end() && i - 2 >= 0 && j - 2 >= 0 && i < seq1.length() && j < seq2.length())
            {
                if (seq2[j] == seq1[i] && seq1[i - 2] == seq2[j - 1] && seq1[i] == seq2[j - 1])
                    return extend_score[make_pair(a, b)];
                else
                    return EXTEND;
            }
            else
                return EXTEND;
        }
        else if (dir == Direction::LEFT)    
        {
            char b = seq1[i - 1];
            char a = '-';
            if (extend_score.find(make_pair(a, b)) != extend_score.end() && i - 2 >= 0 && j - 2 >= 0 && i < seq1.length() && j < seq2.length())
            {
                if (seq2[j] == seq1[i] && seq1[i - 1] == seq2[j - 2] && seq1[i] == seq2[j - 1])
                    return extend_score[make_pair(a, b)];
                else
                    return EXTEND;
            }
            else
                return EXTEND;
        }
    }

    Align(int r, int c, string str1, string str2)
    {
        rows = r;
        cols = c;
        seq1 = str1;
        seq2 = str2;
        middleMat = vector<vector<double>>(rows, vector<double>(cols, INT_MIN));
        xMat = vector<vector<double>>(rows, vector<double>(cols, INT_MIN));
        yMat = vector<vector<double>>(rows, vector<double>(cols, INT_MIN));
        middleDir = vector<vector<char>>(rows, vector<char>(cols, 'm'));
        xDir = vector<vector<char>>(rows, vector<char>(cols, 'x'));
        yDir = vector<vector<char>>(rows, vector<char>(cols, 'y'));
        middleMat[0][0] = 0;
        xMat[0][0] = 0;
        yMat[0][0] = 0;
        Aseq1 = "";
        Aseq2 = "";
        score_m = vector<double>(2);
        score_x = vector<double>(2);
        score_y = vector<double>(2);
    }

private:
    int rows;
    int cols;
    int bandwidth;
    string seq1;
    string seq2;
    vector<vector<double>> middleMat;
    vector<vector<double>> xMat;
    vector<vector<double>> yMat;
    char maxDir;
    vector<vector<char>> middleDir;
    vector<vector<char>> xDir;
    vector<vector<char>> yDir;
    string Aseq1;
    string Aseq2;
    vector<double> score_m;
    vector<double> score_x;
    vector<double> score_y;
};
