#include "fileInputs.cpp"
#include <algorithm>
#include <stack>
#include <unordered_map>
#include <unordered_set>
using namespace std;
class bpRNASequencer
{
public:

    void setDB(string D)
    {
        DB = D;
    }


    void setseq(vector<char> D)
    {
        seq = D;
    }

    void setbp(map <int, int> D)
    {
        bp = D;
    }

    int getNextPair(int idx)
    {
        for (int n = idx + 1; n <= LastPosition; n++)
        {
            if (bp[n] != 0)
            {
                return n;
            }
        }
        return 0;
    }

    int getPrevPair(int idx)
    {
        for (int n = idx - 1; n >= FirstPosition; n--)
        {
            if (bp[n] != 0)
            {
                return n;
            }
        }
        return 0;
    }

    void calculateFirstLastPositions()
    {

        for (auto const& imap : bp)
            keys.push_back(imap.first);
        sort(keys.begin(), keys.end());
        FirstPosition = keys[0];
        LastPosition = keys[keys.size() - 1];

    }

    void getSegments()
    {
        calculateFirstLastPositions();
        int i = getNextPair(0);
        if (i != 0)
        {
            int j = bp[i];
            while (FirstPosition <= i && i <= LastPosition)
            {
                int insegment = 0;
                vector<pair<int, int>> segment;
                if (i < j)
                {

                    segment.push_back(make_pair(i, j));
                    insegment = 1;
                }
                while (insegment)
                {
                    int n = getNextPair(i);
                    int p = getPrevPair(j);
                    if (n && p)
                    {
                        if (bp[n] == p && n < p)
                        {
                            segment.push_back(make_pair(n, p));
                            i = n;
                            j = p;
                        }
                        else
                        {
                            insegment = 0;
                            allsegments.push_back(segment);

                        }
                    }
                }
                i = getNextPair(i);
                j = bp[i];


            }
        }
    }


    int nextPBase(int i)
    {
        while (i < DB.size() && DB[i] == '.')
            i++;
        return i;
    }

    int prevPBase(int i)
    {
        while (i >= 0 && DB[i] == '.')
            i--;
        return i;
    }

    void computePair()
    {
        stack<int> para;
        for (int i = 0; i < DB.size(); i++)
        {
            if (DB[i] == '.')
                bp[i] = -1;
            else if (DB[i] == '(')
            {
                para.push(i);
            }
            else if (DB[i] == ')')
            {
                int p = para.top();
                para.pop();
                bp[i] = p;
                bp[p] = i;
            }
            else
                bp[i] = -1;
        }
    }


    void computeStructureArray()
    {
        char structure;
        for (int i = 0; i < DB.size(); i++)
        {
            if (DB[i] == ')' || DB[i] == '(')
            {
                structure = 'S';
            }
            else
            {
                int nextBPIndex = nextPBase(i + 1);
                int prevBPIndex = prevPBase(i - 1);
                int nextBPPair = bp[nextBPIndex];
                int prevBPPair = bp[prevBPIndex];

                int distance = nextBPIndex - prevBPIndex - 1;
                if (prevBPIndex == -1)
                    structure = 'E';
                else if (DB[prevBPIndex] == '(')
                {
                    if (nextBPIndex == DB.size())
                    {
                        structure = 'E';
                    }
                    else if (DB[nextBPIndex] == '(')
                    {
                        if (nextBPPair == prevBPPair - 1)
                        {
                            structure = 'B';
                        }
                        else
                        {
                            structure = findbtwn(nextBPPair, prevBPPair);
                        }
                    }
                    else if (DB[nextBPIndex] == ')')
                    {
                        structure = 'H';
                    }


                }
                else if (DB[prevBPIndex] == ')')
                {
                    if (nextBPIndex == DB.size())
                    {
                        structure = 'E';
                    }
                    else if (DB[nextBPIndex] == '(')
                        structure = 'X';
                    else if (DB[nextBPIndex] == ')')
                    {

                        if (nextBPPair == prevBPPair - 1)
                        {
                            structure = 'B';
                        }
                        else
                        {
                            structure = findbtwn(nextBPPair, prevBPPair);
                        }

                    }

                }

            }
            structureVec.push_back(structure);
        }
    }


    char findbtwn(int nextBPPair, int prevBPPair)
    {
        nextBPPair++;
        while (nextBPPair < prevBPPair - 1)
        {
            if (DB[nextBPPair] == '(' || DB[nextBPPair] == ')')
            {
                return 'M';
            }
            nextBPPair++;
        }
        return 'I';
    }

    bool pkQuartet(int i, int j, int k, int l)
    {
        if (((i < k) && (k < j) && (j < l)) ||
            ((i < l) && (l < j) && (k < i))) {
            return true;
        }
        return false;
    }

    void separateSegments()
    {
        unordered_map<int, int> knot;
        vector<unordered_set<int>> G(allsegments.size());

        if (allsegments.size() > 1)
        {
            for (size_t i = 0; i < allsegments.size() - 1; i++)
            {
                const vector<pair<int, int>>& segment1 = allsegments[i];
                const pair<int, int>& firstPair1 = segment1.front();
                int s1_5pStart = firstPair1.first;
                int s1_3pStart = firstPair1.second;

                for (size_t j = i + 1; j < allsegments.size(); j++)
                {
                    const vector<pair<int, int>>& segment2 = allsegments[j];
                    const pair<int, int>& firstPair2 = segment2.front();
                    int s2_5pStart = firstPair2.first;
                    int s2_3pStart = firstPair2.second;

                    if (pkQuartet(s1_5pStart, s1_3pStart, s2_5pStart, s2_3pStart))
                    {
                        G[i].insert(j);
                        G[j].insert(i);
                    }
                }
            }
        }

        vector<unordered_set<int>> CC;
        unordered_set<int> visited;

        for (size_t i = 0; i < allsegments.size(); i++)
        {
            if (visited.find(i) == visited.end())
            {
                vector<int> currentCC;
                vector<int> stack;
                stack.push_back(i);

                while (!stack.empty()) {
                    int v = stack.back();
                    stack.pop_back();

                    if (visited.find(v) == visited.end()) {
                        currentCC.push_back(v);
                        visited.insert(v);
                        for (int neighbor : G[v]) {
                            if (visited.find(neighbor) == visited.end()) {
                                stack.push_back(neighbor);
                            }
                        }
                    }
                }

                if (!currentCC.empty())
                {
                    CC.push_back(unordered_set<int>(currentCC.begin(), currentCC.end()));
                }
            }
        }

        vector<vector<pair<int, int>>> segmentsList;
       

        for (const auto& c : CC)
        {
            vector<int> bestKnots;
            string warning;
            bestKnots = getBestKnots(segmentsList, c, G);
            knotsSegment.push_back({});
            for (int v : bestKnots)
            {
                knotsSegment.back().insert(knotsSegment.back().end(), segmentsList[v].begin(), segmentsList[v].end());
                knot[v]++;
            }

        }
        for (int i = 0; i < allsegments.size(); i++)
        {
            if (knot.find(i) != knot.end())
            {
                knotsSegment.push_back(allsegments[i]);
            }
            else
            {
                filteredSegments.push_back(allsegments[i]);
            }
        }

        
    }

    vector<vector<int>> getsubgraph(unordered_set<int> nodes, vector<unordered_set<int>> G) const {
        vector<vector<int>> subG(G.size());
        unordered_set<int> nodesSet(nodes.begin(), nodes.end());

        for (int u : nodesSet) {
            for (int v : nodesSet) {
                if (G[u].count(v) == 1) {
                    subG[u].push_back(v);
                    subG[v].push_back(u);
                }
            }
        }
        return subG;
    }

    vector<int> getBestKnots(vector<vector<pair<int, int>>> segments, unordered_set<int> c, vector<unordered_set<int>> G)
    {
        vector<int> knotsList;
        bool knotleft = 1;
        vector<vector<int>> subGraph = getsubgraph(c, G);
        while (knotleft)
        {
            knotleft = 0;
            vector<vector<int>> CC;
            vector<int> visited;

            for (size_t i = 0; i < segments.size(); i++)
            {
                if (find(visited.begin(), visited.end(), i) != visited.end())
                {
                    vector<int> currentCC;
                    vector<int> stack;
                    stack.push_back(i);

                    while (!stack.empty()) {
                        int v = stack.back();
                        stack.pop_back();

                        if (find(visited.begin(), visited.end(), v) != visited.end()) {
                            currentCC.push_back(v);
                            visited.push_back(v);
                            for (int neighbor : subGraph[v]) {
                                if (find(visited.begin(), visited.end(), neighbor) != visited.end()) {
                                    stack.push_back(neighbor);
                                }
                            }
                        }
                    }

                    if (!currentCC.empty()) {
                        CC.push_back(vector<int>(currentCC.begin(), currentCC.end()));
                    }
                }
            }
            for (const auto& cc : CC)
            {
                int componentSize = cc.size();
                if (componentSize == 2)
                {
                    knotleft = true;

                    auto minV = getMinVPair(cc, segments);
                    subGraph.erase(subGraph.begin(), subGraph.end() + minV);
                    for (int i = 0; i < G.size(); i++)
                    {
                        if (i != minV)
                        {
                            auto itr = find(subGraph[i].begin(), subGraph[i].end(), minV);
                            if (itr != subGraph[i].end())
                            {
                                subGraph[i].erase(itr);

                            }

                        }
                    }
                    knotsList.push_back({ minV });

                }
                else if (componentSize == 3 && isPathGraph(subGraph, cc))
                {

                    int w1 = segments[cc[0]].size();
                    int w2 = segments[cc[1]].size();
                    int w3 = segments[cc[2]].size();

                    if (w2 == w1 + w3) {
                        int v = cc[1] - 1;
                        subGraph.erase(subGraph.begin(), subGraph.end() + v);
                        for (int i = 0; i < G.size(); i++)
                        {
                            if (i != v)
                            {
                                auto itr = find(subGraph[i].begin(), subGraph[i].end(), v);
                                if (itr != subGraph[i].end())
                                {
                                    subGraph[i].erase(itr);
                                }

                            }
                        }
                        knotsList.push_back({ v });
                        knotleft = true;
                    }
                }
                else
                {
                    unordered_map<int, vector<int>> nodeInfo;
                    int minWeight = 10e10;
                    int maxDegree = 0;
                    int maxMaxDegreeScore = -1;
                    int maxDegreeV = -1;
                    int minV = -1;

                    for (int i = 0; i < cc.size(); ++i)
                    {
                        int v = cc[i];
                        int d = G[v].size();
                        int weight = segments[v].size();


                        if (weight < minWeight) {
                            minWeight = weight;
                            minV = v;
                        }


                        if (d >= maxDegree)
                        {
                            if (d > maxDegree) {
                                maxMaxDegreeScore = -1;
                            }

                            int weightSum = 0;
                            vector<int> neighbours(G[v].size());
                            copy(G[v].begin(), G[v].end(), neighbours.begin());

                            for (int w : neighbours) {
                                weightSum += segments[w].size();
                            }

                            int score = weightSum - weight;
                            nodeInfo[v] = { d, score, weight };

                            if (score > maxMaxDegreeScore) {
                                maxDegreeV = v;
                                maxDegree = d;
                                maxMaxDegreeScore = score;
                            }
                        }
                    }


                    int count = 0;
                    for (const auto& entry : nodeInfo) {
                        int v = entry.first;
                        const vector<int>& info = entry.second;
                        int d = info[0];
                        int score = info[1];
                        int weight = info[2];

                        if (d == maxDegree && score == maxMaxDegreeScore) {
                            count++;
                        }
                    }

                    if (maxMaxDegreeScore > 0) {
                        subGraph.erase(subGraph.begin(), subGraph.end() + maxDegreeV);
                        for (int i = 0; i < G.size(); i++)
                        {
                            if (i != maxDegreeV)
                            {
                                auto itr = find(subGraph[i].begin(), subGraph[i].end(), maxDegreeV);
                                if (itr != subGraph[i].end())
                                {
                                    subGraph[i].erase(itr);
                                }

                            }
                        }
                        knotsList.push_back({ maxDegreeV });

                        knotleft = true;



                    }
                    else {
                        subGraph.erase(subGraph.begin(), subGraph.end() + minV);
                        for (int i = 0; i < G.size(); i++)
                        {
                            if (i != minV)
                            {
                                auto itr = find(subGraph[i].begin(), subGraph[i].end(), minV);
                                if (itr != subGraph[i].end())
                                {
                                    subGraph[i].erase(itr);
                                }

                            }
                        }
                        knotsList.push_back({ maxDegreeV });


                        knotleft = true;
                    }
                }
            }




        }
        return knotsList;
    }


    bool isPathGraph(vector<vector<int>> graph, vector<int > cc)
    {
        unordered_set<int> visited;
        visited.insert(cc[0]);

        for (int i = 1; i < cc.size(); ++i)
        {
            if (find(graph[i].begin(), graph[i].end(), cc[i] - 1) != graph[i].end() || visited.find(cc[i]) != visited.end())
            {
                return false;
            }
            visited.insert(cc[i]);
        }

        return true;
    }

    int min(int a, int b)
    {
        return (a < b) ? a : b;
    }

    int getMinVPair(const vector<int>& c, const vector<vector<pair<int, int>>>& segments)
    {
        int v = c[0];
        int w = c[1];
        if (segments[v].size() < segments[w].size())
        {
            return v;
        }
        else if (segments[w].size() < segments[v].size())
        {
            return w;
        }
        else {
            int minV = min(v, w);
            return minV;
        }
    }

    pair<char, char> getBrackets(int n)
    {
        string left = "([{<";
        for (char c = 'A'; c <= 'Z'; c++)
        {
            left += c;
        }

        string right = ")]}>";
        for (char c = 'a'; c <= 'z'; c++)
        {
            right += c;
        }

        return { left[n], right[n] };
    }

    bool knotsCross(const vector<pair<int, int>>& knot1, const vector<pair<int, int>>& knot2)
    {
        vector<pair<int, int>> knot1_copy = knot1;
        int k1_5pStart = knot1_copy[0].first;
        int k1_3pStart = knot1_copy[0].second;

        vector<pair<int, int>> knot2_copy = knot2;
        int k2_5pStart = knot2_copy[0].first;
        int k2_3pStart = knot2_copy[0].second;

        if (pkQuartet(k1_5pStart, k1_3pStart, k2_5pStart, k2_3pStart))
        {
            return true;
        }

        return false;
    }

    void computeDotBracket()
    {
         DB = string(seq.size(), '.');

        for (auto& segment : filteredSegments)
        {
            for (auto& pair : segment)
            {
                int l = pair.first;
                int r = pair.second;
                DB[l - 1] = '(';
                DB[r - 1] = ')';
            }
        }

        for (int i = 0; i < knotsSegment.size(); i++)
        {
            auto brackets = getBrackets(i);

            for (const auto& pair : knotsSegment[i])
            {
                int l = pair.first;
                int r = pair.second;
                DB[l - 1] = brackets.first;
                DB[r - 1] = brackets.second;
            }
        }
    }

    void convertDB()
    {
        computeStructureArray();

    }


    void convertbpSeq()
    {
        getSegments();
        separateSegments();
        computeDotBracket();
        computeStructureArray();
    }

   string getStructureVector()
    {
        return structureVec;
    }

   string getDB()
   {
       return DB;
   }

 private:
    string DB;
    vector<char> seq;
    map <int, int> bp;
    int LastPosition, FirstPosition;
    vector<int> keys;
    vector <vector<pair<int, int>>> allsegments;
    vector<vector<pair<int, int>>> filteredSegments;
    vector<vector<pair<int, int>>> knotsSegment;

    string structureVec;
};


