#include<bits/stdc++.h>

using namespace std;

const long long NMAX = (int)1000 + 5;

int nNodes, nEdges, nDistricts,d[NMAX][NMAX],maxDistance = 0;
vector<pair<int,int>> adjList[NMAX];
double theta = 10.0;
double tau[3] = {0.05,0.05,0.05};
struct BasicUnit{
    int id;
    double x,y;
    double w[3];
}BUs[NMAX];

class Prepare{
public:

    void static readData(){
        ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0);
        cin >> nNodes;
        for(int i = 1;i <= nNodes;i++){
            cin >> BUs[i].id >> BUs[i].x >> BUs[i].y;
            for(int j = 0;j < 3;j++){
                cin >> BUs[i].w[j];
            }
        }
        cin >> nEdges;
        for(int i = 0;i < nEdges;i++){
            int u,v,w;
            cin >> u >> v >> w;
            adjList[u].push_back({v,w});
            adjList[v].push_back({u,w});
        }
        cin >> nDistricts >> tau[0] >> tau[1] >> tau[2];
    }

    void static preProcess(){

        // Compute all-pairs shortest paths using Dijkstra
        // Time complexity: O(n^2 log n + nm log n)

        for(int s = 1;s <= nNodes;s++){
            priority_queue<pair<int,int>,vector<pair<int,int>>,greater<pair<int,int>>> pq;
            for(int i = 1;i <= nNodes;i++){
                d[s][i] = INT_MAX;
            }
            d[s][s] = 0;
            pq.push({0,s});
            while(!pq.empty()){
                auto top = pq.top();
                pq.pop();
                int u = top.second;
                int distU = top.first;
                if(distU > d[s][u]) continue;
                for(auto edge : adjList[u]){
                    int v = edge.first;
                    int w = edge.second;
                    if(d[s][v] > d[s][u] + w){
                        d[s][v] = d[s][u] + w;
                        pq.push({d[s][v],v});
                    }
                }
            }
        }
        for(int i = 1;i <= nNodes;i++){
            for(int j = 1;j <= nNodes;j++){
                if(d[i][j] < INT_MAX && d[i][j] > maxDistance){
                    maxDistance = d[i][j];
                }
            }
        }
    }

    double static meritFunction(vector<int> &assignment,bool verbose = false){
        double compactnessCost = 0.0;
        double balanceCost = 0.0;
        vector<vector<int>> districts(nDistricts + 5,vector<int>());
        double totalWeights[3] = {0,0,0};

        // ✅ Skip unassigned nodes (assignment[i] == -1)
        for(int i = 1;i <= nNodes;i++){
            if(assignment[i] >= 0) {
                districts[assignment[i]].push_back(i);
            }
            totalWeights[0] += BUs[i].w[0];
            totalWeights[1] += BUs[i].w[1];
            totalWeights[2] += BUs[i].w[2];
        }

        totalWeights[0] /= 1.0 * nDistricts;
        totalWeights[1] /= 1.0 * nDistricts;
        totalWeights[2] /= 1.0 * nDistricts;

        for(int i = 0;i < nDistricts;i++){
            // Compute compactness cost: F(X) = (1/dmax) * max diameter
            double districtCompactness = 0.0;
            for(auto u : districts[i]){
                for(auto v : districts[i]){
                    districtCompactness = max(districtCompactness, (double)d[u][v] / maxDistance);
                }
            }
            compactnessCost += districtCompactness;

            // Compute balance cost: G(X) = sum of balance violations
            double currentWeights[3] = {0,0,0};
            for(auto u : districts[i]){
                currentWeights[0] += BUs[u].w[0];
                currentWeights[1] += BUs[u].w[1];
                currentWeights[2] += BUs[u].w[2];
            }
            for(int a = 0; a < 3; a++){
                double lower = (1 - tau[a]) * totalWeights[a];
                double upper = (1 + tau[a]) * totalWeights[a];
                balanceCost += max({0.0, lower - currentWeights[a], currentWeights[a] - upper}) / totalWeights[a];
            }
        }

        if(verbose){
            cerr << setprecision(4) << fixed << "(" << compactnessCost << "," << balanceCost << ")\n";
        }
        return compactnessCost * (1.0) + balanceCost * theta;
    }
};
class ALNS{
private:
    class Operator{
    private: 
        string name;
    public:
        string getName(){
            return name;
        }
        Operator(string _name){
            name = _name;
        }
    };
    class DestroyOperator : public Operator{
    public:
        DestroyOperator(string name) : Operator(name){}
        virtual ~DestroyOperator() = default;
        virtual vector<int> apply(const vector<int> &assignment) = 0;
    };
    class RepairOperator : public Operator{
    public:
        RepairOperator(string name) : Operator(name){}
        virtual ~RepairOperator() = default;
        virtual vector<int> apply(const vector<int> &assignment) = 0;
    };

    static vector<vector<int>> buildDistricts(const vector<int> &assignment){
        vector<vector<int>> districts(nDistricts + 1);
        for(int i = 1; i <= nNodes; i++){
            if(assignment[i] >= 0) districts[assignment[i]].push_back(i);
        }
        return districts;
    }

    static vector<int> getUnassigned(const vector<int> &assignment){
        vector<int> unassigned;
        for(int i = 1; i <= nNodes; i++){
            if(assignment[i] == -1) unassigned.push_back(i);
        }
        return unassigned;
    }

    class RelatedDestroy : public DestroyOperator{
    public:
        RelatedDestroy() : DestroyOperator("RelatedDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            // ✅ Remove nodes near a seed (geographic clustering) - keeps local coherence
            vector<int> res = assignment;
            int seed = 1 + rand() % nNodes;
            vector<pair<int,int>> distList;
            for(int i = 1; i <= nNodes; i++){
                distList.push_back({d[seed][i], i});
            }
            sort(distList.begin(), distList.end());
            int removeCount = max(1, nNodes / 12);
            for(int i = 0; i < removeCount && i < (int)distList.size(); i++){
                res[distList[i].second] = -1;
            }
            return res;
        }
    };

    class CompactnessDestroy : public DestroyOperator{
    public:
        CompactnessDestroy() : DestroyOperator("CompactnessDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            // ✅ Remove nodes that contribute most to diameter in their district
            vector<int> res = assignment;
            vector<pair<double, int>> diameterContrib;
            
            // For each node, compute its max distance to other nodes in same district
            for(int i = 1; i <= nNodes; i++){
                if(res[i] < 0) continue;
                double maxDist = 0;
                for(int j = 1; j <= nNodes; j++){
                    if(res[j] == res[i]){
                        maxDist = max(maxDist, (double)d[i][j]);
                    }
                }
                diameterContrib.push_back({maxDist, i});
            }
            
            sort(diameterContrib.begin(), diameterContrib.end(), greater<pair<double,int>>());
            int removeCount = max(1, nNodes / 15);
            for(int i = 0; i < removeCount && i < (int)diameterContrib.size(); i++){
                res[diameterContrib[i].second] = -1;
            }
            return res;
        }
    };

    class WorstMeritDestroy : public DestroyOperator{
    public:
        WorstMeritDestroy() : DestroyOperator("WorstMeritDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            // ✅ Remove nodes that contribute most to balance violation
            vector<int> res = assignment;
            vector<array<double, 3>> districtWeights(nDistricts, {0.0, 0.0, 0.0});
            
            for(int i = 1; i <= nNodes; i++){
                if(res[i] >= 0){
                    for(int a = 0; a < 3; a++){
                        districtWeights[res[i]][a] += BUs[i].w[a];
                    }
                }
            }
            
            vector<pair<double, int>> impact;
            for(int i = 1; i <= nNodes; i++){
                if(res[i] < 0) continue;
                int k = res[i];
                double violation = 0.0;
                for(int a = 0; a < 3; a++){
                    double newVal = districtWeights[k][a] - BUs[i].w[a];
                    double avg = 0;
                    for(int j = 1; j <= nNodes; j++) avg += BUs[j].w[a];
                    avg = avg / nDistricts;
                    double lower = (1 - tau[a]) * avg;
                    double upper = (1 + tau[a]) * avg;
                    violation += max({0.0, lower - newVal, newVal - upper}) / avg;
                }
                impact.push_back({violation, i});
            }
            
            sort(impact.begin(), impact.end(), greater<pair<double,int>>());
            int removeCount = max(1, nNodes / 15);
            for(int i = 0; i < removeCount && i < (int)impact.size(); i++){
                res[impact[i].second] = -1;
            }
            return res;
        }
    };

    class InterdistrictDestroy : public DestroyOperator{
    public:
        InterdistrictDestroy() : DestroyOperator("InterdistrictDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            // ✅ Remove boundary nodes between districts (encourages district reorganization)
            vector<int> res = assignment;
            vector<pair<int, int>> boundaryNodes;
            
            for(int i = 1; i <= nNodes; i++){
                if(res[i] < 0) continue;
                // Check if node has neighbor in different district
                for(pair<int, int> p : adjList[i]){
                    int j = p.first;
                    if(res[j] >= 0 && res[j] != res[i]){
                        boundaryNodes.push_back({i, res[i]});
                        break;
                    }
                }
            }
            
            int removeCount = max(1, (int)(nNodes * 0.08));
            random_shuffle(boundaryNodes.begin(), boundaryNodes.end());
            for(int i = 0; i < removeCount && i < (int)boundaryNodes.size(); i++){
                res[boundaryNodes[i].first] = -1;
            }
            return res;
        }
    };

    class RandomRepair : public RepairOperator{
    public:
        RandomRepair() : RepairOperator("RandomRepair"){}
        vector<int> apply(const vector<int> &assignment) override{
            // ✅ Assign unassigned to least-loaded district by size
            vector<int> res = assignment;
            vector<int> districtSize(nDistricts, 0);
            for(int i = 1; i <= nNodes; i++){
                if(res[i] >= 0) districtSize[res[i]]++;
            }
            
            for(int i = 1; i <= nNodes; i++){
                if(res[i] == -1){
                    int minDist = min_element(districtSize.begin(), districtSize.end()) - districtSize.begin();
                    res[i] = minDist;
                    districtSize[minDist]++;
                }
            }
            return res;
        }
    };

    class NearestRepair : public RepairOperator{
    public:
        NearestRepair() : RepairOperator("NearestRepair"){}
        vector<int> apply(const vector<int> &assignment) override{
            // ✅ Assign to nearest district (by nearest node in that district)
            // This promotes geographic compactness
            vector<int> res = assignment;
            
            for(int node = 1; node <= nNodes; node++){
                if(res[node] != -1) continue;
                
                int bestDistrict = 0;
                int bestDist = INT_MAX;
                
                for(int k = 0; k < nDistricts; k++){
                    for(int i = 1; i <= nNodes; i++){
                        if(res[i] == k && d[node][i] < bestDist){
                            bestDist = d[node][i];
                            bestDistrict = k;
                        }
                    }
                }
                
                if(bestDist == INT_MAX) bestDistrict = node % nDistricts;
                res[node] = bestDistrict;
            }
            return res;
        }
    };

    class BalancedCompactnessRepair : public RepairOperator{
    public:
        BalancedCompactnessRepair() : RepairOperator("BalancedCompactnessRepair"){}
        vector<int> apply(const vector<int> &assignment) override{
            // ✅ Balance compactness + balance: assign by min(balance_violation + compactness_increase)
            vector<int> res = assignment;
            vector<int> unassigned = getUnassigned(res);
            
            for(int node : unassigned){
                vector<array<double,3>> districtWeights(nDistricts, {0.0, 0.0, 0.0});
                for(int i = 1; i <= nNodes; i++){
                    if(res[i] >= 0){
                        for(int a = 0; a < 3; a++){
                            districtWeights[res[i]][a] += BUs[i].w[a];
                        }
                    }
                }
                
                double bestScore = 1e100;
                int bestDistrict = 0;
                
                for(int k = 0; k < nDistricts; k++){
                    double balanceScore = 0.0;
                    for(int a = 0; a < 3; a++){
                        double newVal = districtWeights[k][a] + BUs[node].w[a];
                        double avg = 0;
                        for(int j = 1; j <= nNodes; j++) avg += BUs[j].w[a];
                        avg = avg / nDistricts;
                        double lower = (1 - tau[a]) * avg;
                        double upper = (1 + tau[a]) * avg;
                        balanceScore += max({0.0, lower - newVal, newVal - upper}) / avg;
                    }
                    
                    // Add compact cost: distance to nearest node in district
                    double compactScore = 1e10;
                    for(int i = 1; i <= nNodes; i++){
                        if(res[i] == k){
                            compactScore = min(compactScore, (double)d[node][i] / maxDistance);
                        }
                    }
                    if(compactScore == 1e10) compactScore = 0; // Empty district
                    
                    double totalScore = balanceScore * 100.0 + compactScore * 1.0;
                    if(totalScore < bestScore){
                        bestScore = totalScore;
                        bestDistrict = k;
                    }
                }
                res[node] = bestDistrict;
            }
            return res;
        }
    };
    static vector<int> greedyInitial(){
        vector<int> assignment(nNodes + 1, -1);
        vector<int> centers;
        set<int> assigned;
        
        // Step 1: Select nDistricts seeds using farthest-first approach
        int seed = 1 + rand() % nNodes;
        centers.push_back(seed);
        assigned.insert(seed);
        assignment[seed] = 0;
        
        // Farthest-first seed selection
        while((int)centers.size() < nDistricts){
            int best = -1;
            int bestDist = -1;
            for(int i = 1; i <= nNodes; i++){
                if(assigned.count(i)) continue;
                int minDist = INT_MAX;
                for(int c : centers){
                    minDist = min(minDist, d[i][c]);
                }
                if(minDist > bestDist){
                    bestDist = minDist;
                    best = i;
                }
            }
            if(best != -1){
                centers.push_back(best);
                assigned.insert(best);
                assignment[best] = centers.size() - 1;
            }
        }
        
        // Step 2: Greedy expansion - assign remaining nodes to nearest center
        for(int i = 1; i <= nNodes; i++){
            if(assigned.count(i)) continue;
            int bestCenter = 0;
            int bestDist = INT_MAX;
            for(int j = 0; j < nDistricts; j++){
                if(d[i][centers[j]] < bestDist){
                    bestDist = d[i][centers[j]];
                    bestCenter = j;
                }
            }
            assignment[i] = bestCenter;
            assigned.insert(i);
        }
        
        return assignment;
    }
public:
    static vector<int> runAlgorithm(int maxIters = 50000, double timeLimit = 3600.0){
        srand((unsigned)time(0));
        
        // Initialize destroy operators - focus on targeted destruction
        vector<shared_ptr<DestroyOperator>> destroyOps = {
            make_shared<RelatedDestroy>(),           // Geographic clustering removal
            make_shared<CompactnessDestroy>(),       // Remove diameter-causing nodes
            make_shared<WorstMeritDestroy>(),        // Remove balance-violation nodes
            make_shared<InterdistrictDestroy>()      // Remove boundary nodes
        };
        
        // Initialize repair operators - focus on targeted construction
        vector<shared_ptr<RepairOperator>> repairOps = {
            make_shared<RandomRepair>(),                    // Distribute to least-loaded
            make_shared<NearestRepair>(),                   // Geographic proximity
            make_shared<BalancedCompactnessRepair>()        // Balance + Compactness trade-off
        };
        
        int nDestroy = destroyOps.size();
        int nRepair = repairOps.size();
        
        // Adaptive weights for operators
        vector<double> destroyWeights(nDestroy, 1.0), repairWeights(nRepair, 1.0);
        vector<double> destroyScores(nDestroy, 0.0), repairScores(nRepair, 0.0);
        vector<double> destroyCounts(nDestroy, 1e-6), repairCounts(nRepair, 1e-6);
        vector<double> scores = {0.5, 2.0, 10.0};  // small/medium/best reward
        
        // Initialize current and best solution using greedy algorithm
        vector<int> curAssignment = greedyInitial();
        vector<int> bestAssignment = curAssignment;
        
        double curMerit = Prepare::meritFunction(curAssignment);
        double bestMerit = curMerit;
        
        // Simulated Annealing parameters
        double T = max(1.0, curMerit);
        double alpha = 0.9995;
        
        auto startTime = chrono::high_resolution_clock::now();
        int noImprove = 0;
        
        // Main ALNS loop
        for(int iteration = 0; iteration < maxIters; iteration++){
            
            // Check time limit
            auto elapsed = chrono::duration<double>(chrono::high_resolution_clock::now() - startTime).count();
            if(elapsed > timeLimit) break;
            
            // Roulette selection of operators
            int destroyIdx = rouletteSelect(destroyWeights);
            int repairIdx = rouletteSelect(repairWeights);
            // Perform destroy
            vector<int> destroyed = destroyOps[destroyIdx]->apply(curAssignment);
            
            // Perform repair
            vector<int> repaired = repairOps[repairIdx]->apply(destroyed);
            
            // Evaluate new solution
            cerr << "\n#Iter " << iteration << ": Destroy(" << destroyOps[destroyIdx]->getName() << ") Repair(" << repairOps[repairIdx]->getName() << "): \n";
            double newMerit = Prepare::meritFunction(repaired,true);
            cerr << "Current Merit: " << curMerit << ", New Merit: " << newMerit << ", Best Merit: " << bestMerit << "\n";
            // cerr << "Iteration " << iteration << " Repaired Merit: ";
            // cerr << Prepare::meritFunction(repaired, true) << "\n";
            // cerr << "Iter " << iteration << ": Destroy(" << destroyOps[destroyIdx]->getName() << ") Repair(" << repairOps[repairIdx]->getName() << ") New Merit: " << newMerit << " Current Merit: " << curMerit << " Best Merit: " << bestMerit << "\n";
            // SA acceptance criterion
            bool accepted = false;
            if(newMerit <= curMerit){
                accepted = true;
            } else {
                double prob = exp((curMerit - newMerit) / max(1e-9, T));
                if((double)rand() / RAND_MAX < prob) accepted = true;
            }
            
            if(accepted){
                curAssignment = repaired;
                curMerit = newMerit;
                
                // Check if best solution improved (lexicographic or merit)
                if(newMerit < bestMerit){
                    bestAssignment = repaired;
                    bestMerit = newMerit;
                    destroyScores[destroyIdx] += scores[2];
                    repairScores[repairIdx] += scores[2];
                    noImprove = 0;
                } else {
                    destroyScores[destroyIdx] += scores[1];
                    repairScores[repairIdx] += scores[1];
                }
            } else {
                destroyScores[destroyIdx] += scores[0];
                repairScores[repairIdx] += scores[0];
            }
            
            destroyCounts[destroyIdx]++;
            repairCounts[repairIdx]++;
            
            // Periodically update adaptive weights
            if(iteration % 50 == 0 && iteration > 0){
                for(int i = 0; i < nDestroy; i++){
                    destroyWeights[i] = 0.8 * destroyWeights[i] + 0.2 * (destroyScores[i] / destroyCounts[i]);
                }
                for(int i = 0; i < nRepair; i++){
                    repairWeights[i] = 0.8 * repairWeights[i] + 0.2 * (repairScores[i] / repairCounts[i]);
                }
            }
            
            // Cool down temperature
            T *= alpha;
            noImprove++;
            
            // Restart if stagnant
            if(noImprove > 500){
                curAssignment = greedyInitial();
                curMerit = Prepare::meritFunction(curAssignment);
                noImprove = 0;
            }
        }
        
        return bestAssignment;
    }
    
    static int rouletteSelect(const vector<double> &weights){
        double sum = 0;
        for(auto w : weights) sum += w;
        if(sum == 0) return rand() % weights.size();
        double r = ((double)rand() / RAND_MAX) * sum;
        double cum = 0;
        for(int i = 0; i < (int)weights.size(); i++){
            cum += weights[i];
            if(r <= cum) return i;
        }
        return weights.size() - 1;
    }
};
signed main(){

    if(ifstream("input.txt")){
        freopen("input.txt","r",stdin);
        freopen("output.txt","w",stdout);
    }

    Prepare::readData();
    Prepare::preProcess();
    
    cerr << "Running ALNS Algorithm...\n";
    vector<int> bestAssignment = ALNS::runAlgorithm(50000, 3600.0);
    
    // Output results
    double bestMerit = Prepare::meritFunction(bestAssignment);
    cerr << Prepare::meritFunction(bestAssignment, true) << "\n";
    cerr << "Best Merit: " << bestMerit << "\n";
    cerr << "Solution: ";
    for(int i = 1; i <= nNodes; i++){
        cerr << bestAssignment[i] << " ";
        cout << bestAssignment[i] << " ";
    }
    cerr << "\n";
    cout << "\n";

    
    return 0;
}