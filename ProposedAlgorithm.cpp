#include<bits/stdc++.h>

using namespace std;

mt19937 rng((unsigned)chrono::steady_clock::now().time_since_epoch().count());

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
        for(int i = 0;i < nNodes;i++){
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

        for(int s = 0;s < nNodes;s++){
            priority_queue<pair<int,int>,vector<pair<int,int>>,greater<pair<int,int>>> pq;
            for(int i = 0;i < nNodes;i++){
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
                    long long newDist = (long long)d[s][u] + w;
                    if(newDist < d[s][v]){
                        d[s][v] = (int)newDist;
                        pq.push({d[s][v],v});
                    }
                }
            }
        }
        for(int i = 0;i < nNodes;i++){
            for(int j = 0;j < nNodes;j++){
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
        for(int i = 0;i < nNodes;i++){
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
                    if(maxDistance > 0)
                        districtCompactness = max(districtCompactness, (double)d[u][v] / maxDistance);
                }
            }
            compactnessCost = max(compactnessCost, districtCompactness);

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
                if(totalWeights[a] > 1e-12)
                    balanceCost += max({0.0, lower - currentWeights[a], currentWeights[a] - upper}) / totalWeights[a];
            }
        }

        if(verbose){
            cerr << setprecision(4) << fixed << "(" << compactnessCost << "," << balanceCost << ")\n";
        }
        return compactnessCost * (1.0) + balanceCost * theta;
    }
};

// ====================== LP SOLVER (Two-Phase Simplex) ======================
// Giải: maximize c^T x, s.t. Ax <= b, x >= 0
struct LPSolver {
    int m, n;
    vector<int> N, B;
    vector<vector<double>> D;

    LPSolver(vector<vector<double>>& A, vector<double>& b, vector<double>& c) :
        m(A.size()), n(A[0].size()), N(n + 1), B(m),
        D(m + 2, vector<double>(n + 2, 0.0)) {
        for(int i = 0; i < m; i++)
            for(int j = 0; j < n; j++)
                D[i][j] = A[i][j];
        for(int i = 0; i < m; i++){
            B[i] = n + i;
            D[i][n] = -1;
            D[i][n + 1] = b[i];
        }
        for(int j = 0; j < n; j++){
            N[j] = j;
            D[m][j] = -c[j];
        }
        N[n] = -1;
        D[m + 1][n] = 1;
    }

    void pivot(int r, int s){
        double inv = 1.0 / D[r][s];
        for(int i = 0; i <= m + 1; i++) if(i != r){
            for(int j = 0; j <= n + 1; j++) if(j != s)
                D[i][j] -= D[i][s] * D[r][j] * inv;
            D[i][s] *= -inv;
        }
        for(int j = 0; j <= n + 1; j++) if(j != s) D[r][j] *= inv;
        D[r][s] = inv;
        swap(B[r], N[s]);
    }

    bool simplex(int phase){
        int x = phase == 1 ? m + 1 : m;
        while(true){
            int s = -1;
            for(int j = 0; j <= n; j++){
                if(phase == 2 && N[j] == -1) continue;
                if(s == -1 || D[x][j] < D[x][s] || (D[x][j] == D[x][s] && N[j] < N[s]))
                    s = j;
            }
            if(D[x][s] > -1e-9) return true;
            int r = -1;
            for(int i = 0; i < m; i++){
                if(D[i][s] < 1e-9) continue;
                if(r == -1 || D[i][n+1] / D[i][s] < D[r][n+1] / D[r][s] ||
                   (D[i][n+1] / D[i][s] == D[r][n+1] / D[r][s] && B[i] < B[r]))
                    r = i;
            }
            if(r == -1) return false; // unbounded
            pivot(r, s);
        }
    }

    double solve(vector<double>& x){
        int r = 0;
        for(int i = 1; i < m; i++)
            if(D[i][n+1] < D[r][n+1]) r = i;
        if(D[r][n+1] < -1e-9){
            pivot(r, n);
            if(!simplex(1) || D[m+1][n+1] < -1e-9)
                return -numeric_limits<double>::infinity();
            for(int i = 0; i < m; i++)
                if(B[i] == -1){
                    int s = -1;
                    for(int j = 0; j <= n; j++)
                        if(s == -1 || D[i][j] < D[i][s] || (D[i][j] == D[i][s] && N[j] < N[s]))
                            s = j;
                    pivot(i, s);
                }
        }
        bool feas = simplex(2);
        x.assign(n, 0);
        for(int i = 0; i < m; i++)
            if(B[i] < n) x[B[i]] = D[i][n+1];
        return feas ? D[m][n+1] : -numeric_limits<double>::infinity();
    }
};

class ALNS{
private:
    class Operator{
    private: 
        string name;
    public:
        string getName(){ return name; }
        Operator(string _name) : name(_name){}
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
        vector<vector<int>> districts(nDistricts);
        for(int i = 0; i < nNodes; i++){
            if(assignment[i] >= 0) districts[assignment[i]].push_back(i);
        }
        return districts;
    }

    static vector<int> getUnassigned(const vector<int> &assignment){
        vector<int> unassigned;
        for(int i = 0; i < nNodes; i++){
            if(assignment[i] == -1) unassigned.push_back(i);
        }
        return unassigned;
    }

    // Tính medoid (node có tổng khoảng cách đến các node khác trong cùng district nhỏ nhất)
    static int findMedoid(const vector<int> &members){
        if(members.empty()) return -1;
        int best = members[0];
        long long bestSum = LLONG_MAX;
        for(int u : members){
            long long s = 0;
            for(int v : members) s += d[u][v];
            if(s < bestSum){ bestSum = s; best = u; }
        }
        return best;
    }

    // Tính diameter của 1 district
    static int districtDiameter(const vector<int> &members){
        int diam = 0;
        for(int i = 0; i < (int)members.size(); i++)
            for(int j = i+1; j < (int)members.size(); j++)
                diam = max(diam, d[members[i]][members[j]]);
        return diam;
    }

    // Tính max distance từ node u đến bất kỳ node nào trong district
    static int maxDistTo(int u, const vector<int> &members){
        int mx = 0;
        for(int v : members) mx = max(mx, d[u][v]);
        return mx;
    }

    // ====================== DESTROY OPERATORS ======================

    // Destroy 1: DiameterDestroy
    // Tìm district có diameter lớn nhất → xóa các node xa medoid nhất
    // → trực tiếp giảm objective (max diameter)
    class DiameterDestroy : public DestroyOperator{
    public:
        DiameterDestroy() : DestroyOperator("DiameterDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            auto districts = buildDistricts(res);
            // tìm district có diameter lớn nhất
            int worstK = 0, worstDiam = 0;
            for(int k = 0; k < nDistricts; k++){
                int diam = districtDiameter(districts[k]);
                if(diam > worstDiam){ worstDiam = diam; worstK = k; }
            }
            // tìm medoid của district đó
            int med = findMedoid(districts[worstK]);
            // sắp xếp node theo khoảng cách đến medoid giảm dần
            vector<pair<int,int>> distToMed;
            for(int u : districts[worstK])
                distToMed.push_back({d[u][med], u});
            sort(distToMed.rbegin(), distToMed.rend());
            // xóa 20-30% node xa nhất (nhưng giữ ít nhất 1 node)
            int removeCount = max(1, (int)(districts[worstK].size() * 0.25));
            removeCount = min(removeCount, (int)distToMed.size() - 1);
            for(int i = 0; i < removeCount; i++)
                res[distToMed[i].second] = -1;
            return res;
        }
    };

    // Destroy 2: WorstBalanceDestroy
    // Xóa node gây balance violation lớn nhất
    class WorstBalanceDestroy : public DestroyOperator{
    public:
        WorstBalanceDestroy() : DestroyOperator("WorstBalanceDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            double avg[3] = {0,0,0};
            for(int i = 0; i < nNodes; i++)
                for(int a = 0; a < 3; a++) avg[a] += BUs[i].w[a];
            for(int a = 0; a < 3; a++) avg[a] /= nDistricts;

            vector<array<double,3>> dw(nDistricts, {0,0,0});
            for(int i = 0; i < nNodes; i++)
                if(res[i] >= 0)
                    for(int a = 0; a < 3; a++) dw[res[i]][a] += BUs[i].w[a];

            // Tìm district vi phạm balance nhiều nhất
            int worstK = 0;
            double worstViol = 0;
            for(int k = 0; k < nDistricts; k++){
                double v = 0;
                for(int a = 0; a < 3; a++){
                    if(avg[a] < 1e-12) continue;
                    double lo = (1 - tau[a]) * avg[a], hi = (1 + tau[a]) * avg[a];
                    v += max({0.0, lo - dw[k][a], dw[k][a] - hi}) / avg[a];
                }
                if(v > worstViol){ worstViol = v; worstK = k; }
            }
            // Xóa node trong district đó gây tăng violation many nhất
            auto districts = buildDistricts(res);
            vector<pair<double,int>> impact;
            for(int u : districts[worstK]){
                double costWithout = 0;
                for(int a = 0; a < 3; a++){
                    double nw = dw[worstK][a] - BUs[u].w[a];
                    double lo = (1 - tau[a]) * avg[a], hi = (1 + tau[a]) * avg[a];
                    costWithout += max({0.0, lo - nw, nw - hi}) / avg[a];
                }
                // improvement khi xóa u
                impact.push_back({worstViol - costWithout, u});
            }
            sort(impact.rbegin(), impact.rend());
            int removeCount = max(1, nNodes / 15);
            removeCount = min(removeCount, (int)impact.size() - 1);
            for(int i = 0; i < removeCount; i++)
                res[impact[i].second] = -1;
            return res;
        }
    };

    // Destroy 3: RandomDestroy
    // Xóa ngẫu nhiên 5-10% node → đa dạng hóa tìm kiếm
    class RandomDestroy : public DestroyOperator{
    public:
        RandomDestroy() : DestroyOperator("RandomDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            vector<int> nodes;
            for(int i = 0; i < nNodes; i++) nodes.push_back(i);
            shuffle(nodes.begin(), nodes.end(), rng);
            int removeCount = max(1, (int)(nNodes * 0.07));
            for(int i = 0; i < removeCount; i++)
                res[nodes[i]] = -1;
            return res;
        }
    };

    // Destroy 4: DistrictGroupDestroy
    // Chọn 2-3 district có diameter xấu nhất, xóa toàn bộ trừ medoid
    class DistrictGroupDestroy : public DestroyOperator{
    public:
        DistrictGroupDestroy() : DestroyOperator("DistrictGroupDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            auto districts = buildDistricts(res);
            // tính diameter từng district
            vector<pair<int,int>> diamList; // (diameter, district_id)
            for(int k = 0; k < nDistricts; k++)
                diamList.push_back({districtDiameter(districts[k]), k});
            sort(diamList.rbegin(), diamList.rend());
            // chọn 2-3 district có diameter lớn nhất
            int groupSize = 2 + rng() % 2; // 2 hoặc 3
            groupSize = min(groupSize, nDistricts);
            for(int i = 0; i < groupSize; i++){
                int k = diamList[i].second;
                int med = findMedoid(districts[k]);
                for(int u : districts[k])
                    if(u != med) res[u] = -1;
            }
            return res;
        }
    };

    // Destroy 5: DiameterPairDestroy
    // Tìm cặp node xa nhất trong cùng district, xóa cả 2 + các node lân cận chúng
    class DiameterPairDestroy : public DestroyOperator{
    public:
        DiameterPairDestroy() : DestroyOperator("DiameterPairDestroy"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            auto districts = buildDistricts(res);
            // tìm cặp xa nhất toàn cục
            int worstU = -1, worstV = -1, worstDist = 0;
            for(int k = 0; k < nDistricts; k++){
                for(int i = 0; i < (int)districts[k].size(); i++)
                    for(int j = i+1; j < (int)districts[k].size(); j++){
                        int dd = d[districts[k][i]][districts[k][j]];
                        if(dd > worstDist){
                            worstDist = dd;
                            worstU = districts[k][i];
                            worstV = districts[k][j];
                        }
                    }
            }
            if(worstU == -1) return res;
            // xóa worstU, worstV, + vài node gần chúng
            res[worstU] = -1;
            res[worstV] = -1;
            // xóa thêm node xa center trong cùng district
            int dk = assignment[worstU];
            int med = findMedoid(districts[dk]);
            vector<pair<int,int>> dists;
            for(int u : districts[dk])
                if(u != med) dists.push_back({d[u][med], u});
            sort(dists.rbegin(), dists.rend());
            int extra = max(1, (int)dists.size() / 5);
            for(int i = 0; i < extra && i < (int)dists.size(); i++)
                res[dists[i].second] = -1;
            return res;
        }
    };

    // ====================== REPAIR OPERATORS ======================

    // Hàm tính insertion cost: balance_penalty * BIG + diameter_increase
    // Không cần liên thông → xét TẤT CẢ district
    static double insertionCost(int node, int k, 
                                const vector<int> &res,
                                const vector<array<double,3>> &dw,
                                const double avg[3],
                                const vector<vector<int>> &districts){
        // balance violation delta
        double balDelta = 0;
        for(int a = 0; a < 3; a++){
            if(avg[a] < 1e-12) continue;
            double oldW = dw[k][a];
            double newW = oldW + BUs[node].w[a];
            double lo = (1 - tau[a]) * avg[a], hi = (1 + tau[a]) * avg[a];
            double oldV = max({0.0, lo - oldW, oldW - hi});
            double newV = max({0.0, lo - newW, newW - hi});
            balDelta += (newV - oldV) / avg[a];
        }
        // diameter increase: max distance từ node đến bất kỳ node nào trong district k
        double diamInc = 0;
        if(maxDistance > 0)
            for(int v : districts[k])
                diamInc = max(diamInc, (double)d[node][v] / maxDistance);
        return balDelta * 1000.0 + diamInc;
    }

    // Repair 1: NearestCenterRepair
    // Gán mỗi node chưa gán vào district có medoid gần nhất
    class NearestCenterRepair : public RepairOperator{
    public:
        NearestCenterRepair() : RepairOperator("NearestCenterRepair"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            auto districts = buildDistricts(res);
            // tìm medoid mỗi district
            vector<int> medoids(nDistricts);
            for(int k = 0; k < nDistricts; k++)
                medoids[k] = findMedoid(districts[k]);
            // gán
            for(int node = 0; node < nNodes; node++){
                if(res[node] != -1) continue;
                int bestK = 0, bestD = INT_MAX;
                for(int k = 0; k < nDistricts; k++){
                    if(medoids[k] == -1) continue;
                    if(d[node][medoids[k]] < bestD){
                        bestD = d[node][medoids[k]];
                        bestK = k;
                    }
                }
                res[node] = bestK;
                districts[bestK].push_back(node);
            }
            return res;
        }
    };

    // Repair 2: MinDiameterRepair
    // Gán mỗi node vào district mà max distance đến node hiện tại là nhỏ nhất
    // → trực tiếp minimize diameter increase
    class MinDiameterRepair : public RepairOperator{
    public:
        MinDiameterRepair() : RepairOperator("MinDiameterRepair"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            auto districts = buildDistricts(res);
            vector<int> unassigned = getUnassigned(res);
            // sắp xếp unassigned: ưu tiên node "khó" (xa mọi center nhất)
            vector<int> medoids(nDistricts);
            for(int k = 0; k < nDistricts; k++)
                medoids[k] = findMedoid(districts[k]);
            sort(unassigned.begin(), unassigned.end(), [&](int a, int b){
                int minA = INT_MAX, minB = INT_MAX;
                for(int k = 0; k < nDistricts; k++){
                    if(medoids[k] >= 0){ minA = min(minA, d[a][medoids[k]]); minB = min(minB, d[b][medoids[k]]); }
                }
                return minA > minB; // node xa nhất xếp trước
            });
            for(int node : unassigned){
                int bestK = 0;
                int bestMaxDist = INT_MAX;
                for(int k = 0; k < nDistricts; k++){
                    if(districts[k].empty()) continue; // Bỏ qua district rỗng
                    int mx = maxDistTo(node, districts[k]);
                    if(mx < bestMaxDist){ bestMaxDist = mx; bestK = k; }
                }
                res[node] = bestK;
                districts[bestK].push_back(node);
            }
            return res;
        }
    };

    // Repair 3: BalancedCompactnessRepair
    // Gán theo tradeoff balance + compactness. Không yêu cầu liên thông
    class BalancedCompactnessRepair : public RepairOperator{
    public:
        BalancedCompactnessRepair() : RepairOperator("BalancedCompactnessRepair"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            auto districts = buildDistricts(res);
            vector<int> unassigned = getUnassigned(res);

            vector<array<double,3>> dw(nDistricts, {0,0,0});
            double avg[3] = {0,0,0};
            for(int i = 0; i < nNodes; i++)
                for(int a = 0; a < 3; a++) avg[a] += BUs[i].w[a];
            for(int a = 0; a < 3; a++) avg[a] /= nDistricts;
            for(int i = 0; i < nNodes; i++)
                if(res[i] >= 0) for(int a = 0; a < 3; a++) dw[res[i]][a] += BUs[i].w[a];

            for(int node : unassigned){
                double bestScore = 1e100;
                int bestK = 0;
                for(int k = 0; k < nDistricts; k++){
                    double cost = insertionCost(node, k, res, dw, avg, districts);
                    if(cost < bestScore){ bestScore = cost; bestK = k; }
                }
                res[node] = bestK;
                districts[bestK].push_back(node);
                for(int a = 0; a < 3; a++) dw[bestK][a] += BUs[node].w[a];
            }
            return res;
        }
    };

    // Repair 4: RegretRepair (2-regret)
    // Gán node có regret lớn nhất (= chênh lệch cost giữa choice tốt nhất và nhì)
    // Ưu tiên node "quan trọng" — ít lựa chọn tốt → giảm rủi ro gán sai
    class RegretRepair : public RepairOperator{
    public:
        RegretRepair() : RepairOperator("RegretRepair"){}
        vector<int> apply(const vector<int> &assignment) override{
            vector<int> res = assignment;
            auto districts = buildDistricts(res);
            vector<int> unassigned = getUnassigned(res);

            vector<array<double,3>> dw(nDistricts, {0,0,0});
            double avg[3] = {0,0,0};
            for(int i = 0; i < nNodes; i++)
                for(int a = 0; a < 3; a++) avg[a] += BUs[i].w[a];
            for(int a = 0; a < 3; a++) avg[a] /= nDistricts;
            for(int i = 0; i < nNodes; i++)
                if(res[i] >= 0) for(int a = 0; a < 3; a++) dw[res[i]][a] += BUs[i].w[a];

            while(!unassigned.empty()){
                double bestRegret = -1e100;
                int bestNode = -1, bestK = 0;
                double bestFirst = 1e100;

                for(int node : unassigned){
                    vector<pair<double,int>> costs;
                    for(int k = 0; k < nDistricts; k++)
                        costs.push_back({insertionCost(node, k, res, dw, avg, districts), k});
                    sort(costs.begin(), costs.end());
                    double regret = (costs.size() > 1 ? costs[1].first : 1e15) - costs[0].first;
                    if(regret > bestRegret || (regret == bestRegret && costs[0].first < bestFirst)){
                        bestRegret = regret;
                        bestNode = node;
                        bestK = costs[0].second;
                        bestFirst = costs[0].first;
                    }
                }
                if(bestNode == -1) break;
                res[bestNode] = bestK;
                districts[bestK].push_back(bestNode);
                for(int a = 0; a < 3; a++) dw[bestK][a] += BUs[bestNode].w[a];
                unassigned.erase(find(unassigned.begin(), unassigned.end(), bestNode));
            }
            return res;
        }
    };

    // ====================== INITIAL SOLUTION ======================

    // Tham lam: chọn p seed ngẫu nhiên (mỗi district 1 node),
    // sau đó thêm từng node vào district mà merit tăng ít nhất.
    static vector<int> greedyInitial(){
        vector<int> assignment(nNodes, -1);
        
        // Chọn p seed ngẫu nhiên
        vector<int> allNodes;
        for(int i = 0; i < nNodes; i++) allNodes.push_back(i);
        shuffle(allNodes.begin(), allNodes.end(), rng);
        
        for(int k = 0; k < nDistricts && k < nNodes; k++)
            assignment[allNodes[k]] = k;
        
        // Danh sách node chưa gán
        vector<int> unassigned;
        for(int i = 0; i < nNodes; i++)
            if(assignment[i] == -1) unassigned.push_back(i);
        
        // Shuffle để tránh bias thứ tự
        shuffle(unassigned.begin(), unassigned.end(), rng);
        
        // Thêm từng node vào district mà merit tăng ít nhất
        for(int node : unassigned){
            double bestMerit = 1e100;
            int bestK = 0;
            for(int k = 0; k < nDistricts; k++){
                assignment[node] = k;
                double m = Prepare::meritFunction(assignment);
                if(m < bestMerit){ bestMerit = m; bestK = k; }
            }
            assignment[node] = bestK;
        }
        
        return assignment;
    }

    // ====================== SIMPLEX-BASED INITIAL SOLUTION ======================
    // Giải LP relaxation của mô hình DTDP bằng Simplex đơn hình hai pha
    //
    // min  D
    // s.t. D >= d_ij * (x_im + x_jm - 1)    ∀i,j ∈ V, m ∈ {1..p}
    //      Σ_m x_im = 1                      ∀i ∈ V
    //      (1/μ^a) Σ_i w_i^a x_im >= (1-τ^a) ∀a, m
    //      (1/μ^a) Σ_i w_i^a x_im <= (1+τ^a) ∀a, m
    //      x_im ∈ [0,1]
    //
    // Relax: bỏ ràng buộc diameter, chỉ giữ gán + cân bằng
    // → LP nhỏ, giải nhanh, tạo lời giải feasible về balance
    // Sau đó round fractional solution → integer assignment

    static vector<int> simplexInitial(){
        int V = nNodes;
        int P = nDistricts;
        int nVars = V * P; // x_im: (i-1)*P + m

        // Giới hạn kích thước LP — quá lớn thì fallback greedy
        if((long long)V * P > 8000){
            cerr << "[Simplex] LP qua lon (" << V << "x" << P << "=" << V*P 
                 << " vars), fallback greedy.\n";
            return greedyInitial();
        }

        cerr << "[Simplex] Khoi tao LP: " << nVars << " bien (chi balance)...\n";

        // idx(i, m) = i*P + m, i ∈ [0..V-1], m ∈ [0..P-1]
        auto idx = [&](int i, int m) -> int { return i * P + m; };

        // Tính μ^a = (Σ w_i^a) / P
        double mu[3] = {0, 0, 0};
        for(int i = 0; i < V; i++)
            for(int a = 0; a < 3; a++)
                mu[a] += BUs[i].w[a];
        for(int a = 0; a < 3; a++) mu[a] /= P;

        // ---- Xây dựng ràng buộc: gán + cân bằng ----
        vector<vector<double>> A;
        vector<double> b;

        // (1) Σ_m x_im <= 1, ∀i
        for(int i = 0; i < V; i++){
            vector<double> row(nVars, 0.0);
            for(int m = 0; m < P; m++) row[idx(i, m)] = 1.0;
            A.push_back(row);
            b.push_back(1.0);
        }
        // (2) -Σ_m x_im <= -1 (tức Σ_m x_im >= 1), ∀i
        for(int i = 0; i < V; i++){
            vector<double> row(nVars, 0.0);
            for(int m = 0; m < P; m++) row[idx(i, m)] = -1.0;
            A.push_back(row);
            b.push_back(-1.0);
        }
        // (3) Balance upper: Σ_i w_i^a x_im <= μ^a(1+τ^a), ∀a ∈ {0,1,2}, m ∈ [0..P-1]
        for(int a = 0; a < 3; a++){
            for(int m = 0; m < P; m++){
                vector<double> row(nVars, 0.0);
                for(int i = 0; i < V; i++)
                    row[idx(i, m)] = BUs[i].w[a];
                A.push_back(row);
                b.push_back(mu[a] * (1 + tau[a]));
            }
        }
        // (4) Balance lower: -Σ_i w_i^a x_im <= -μ^a(1-τ^a), ∀a, m
        for(int a = 0; a < 3; a++){
            for(int m = 0; m < P; m++){
                vector<double> row(nVars, 0.0);
                for(int i = 0; i < V; i++)
                    row[idx(i, m)] = -BUs[i].w[a];
                A.push_back(row);
                b.push_back(-mu[a] * (1 - tau[a]));
            }
        }

        // Objective: maximize 0 (feasibility LP — chỉ cần tìm nghiệm feasible)
        // Dùng heuristic: maximize Σ x_im / dist(i, center_m)
        // để kết quả compact hơn khi round
        vector<double> c(nVars, 0.0);

        // Giải LP
        cerr << "[Simplex] Constraints = " << A.size() << "\n";
        LPSolver lp(A, b, c);
        vector<double> xopt;
        double obj = lp.solve(xopt);

        if(obj == -numeric_limits<double>::infinity() || xopt.empty()){
            cerr << "[Simplex] LP infeasible, fallback greedy.\n";
            return greedyInitial();
        }

        // ---- Round fractional solution → integer assignment ----
        // Gán node i vào district m có x_im lớn nhất
        vector<int> assignment(V, -1);
        for(int i = 0; i < V; i++){
            int bestM = 0;
            double bestVal = -1;
            for(int m = 0; m < P; m++){
                int idxVal = idx(i, m);
                double val = (idxVal < (int)xopt.size()) ? xopt[idxVal] : 0.0;
                if(val > bestVal){ bestVal = val; bestM = m; }
            }
            assignment[i] = bestM;
        }

        cerr << "[Simplex] Done. Feasible balance solution found.\n";
        return assignment;
    }
public:
    static vector<int> runAlgorithm(int maxIters = 50000, double timeLimit = 3600.0){
        // rng đã khởi tạo ở global scope
        
        vector<shared_ptr<DestroyOperator>> destroyOps = {
            make_shared<DiameterDestroy>(),          // Xóa node xa medoid trong district có diameter lớn nhất
            make_shared<WorstBalanceDestroy>(),       // Xóa node gây balance violation
            make_shared<RandomDestroy>(),             // Xóa ngẫu nhiên → đa dạng hóa
            make_shared<DistrictGroupDestroy>(),      // Reset 2-3 district tệ nhất
            make_shared<DiameterPairDestroy>()        // Xóa cặp node tạo diameter + lân cận
        };
        
        vector<shared_ptr<RepairOperator>> repairOps = {
            make_shared<NearestCenterRepair>(),       // Gán vào district gần medoid nhất
            make_shared<MinDiameterRepair>(),         // Minimize max distance khi gán
            make_shared<BalancedCompactnessRepair>(), // Balance + compactness trade-off
            make_shared<RegretRepair>()               // 2-Regret: ưu tiên node khó
        };
        
        int nDestroy = destroyOps.size();
        int nRepair = repairOps.size();
        
        vector<double> destroyWeights(nDestroy, 1.0), repairWeights(nRepair, 1.0);
        vector<double> destroyScores(nDestroy, 0.0), repairScores(nRepair, 0.0);
        vector<double> destroyCounts(nDestroy, 1e-6), repairCounts(nRepair, 1e-6);
        vector<double> scores = {0.5, 2.0, 10.0};
        
        // Giải simplex 1 lần, lưu kết quả để tái sử dụng khi restart
        vector<int> cachedSimplexSol = simplexInitial();
        double simplexMerit = Prepare::meritFunction(cachedSimplexSol);
        cerr << "[Init] Simplex merit = " << simplexMerit << "\n";

        vector<int> greedySol = greedyInitial();
        cerr << "Assignement after greedy initialization:\n";
        for(int i = 0; i < nNodes; i++){
            cerr << greedySol[i] << " ";
        }
        double greedyMerit = Prepare::meritFunction(greedySol);
        cerr << "[Init] Greedy merit = " << greedyMerit << "\n";

        vector<int> curAssignment = (greedyMerit < simplexMerit) ? greedySol : cachedSimplexSol;
        vector<int> bestAssignment = curAssignment;
        //cerr << "[Init] Chon: Greedy\n";
        
        double curMerit = Prepare::meritFunction(curAssignment);
        double bestMerit = curMerit;
        
        double T = max(1.0, curMerit * 0.05);
        double alpha = 0.9995;
        
        auto startTime = chrono::high_resolution_clock::now();
        int noImprove = 0;
        int restartCount = 0;
        
        for(int iteration = 0; iteration < maxIters; iteration++){
            auto elapsed = chrono::duration<double>(chrono::high_resolution_clock::now() - startTime).count();
            if(elapsed > timeLimit) break;
            
            int destroyIdx = rouletteSelect(destroyWeights);
            int repairIdx = rouletteSelect(repairWeights);
            
            vector<int> destroyed = destroyOps[destroyIdx]->apply(curAssignment);
            vector<int> repaired = repairOps[repairIdx]->apply(destroyed);
            
            double newMerit = Prepare::meritFunction(repaired);

            bool accepted = false;
            if(newMerit <= curMerit){
                accepted = true;
            } else {
                double prob = exp((curMerit - newMerit) / max(1e-9, T));
                uniform_real_distribution<double> dist01(0.0, 1.0);
                if(dist01(rng) < prob) accepted = true;
            }
            
            if(accepted){
                curAssignment = repaired;
                curMerit = newMerit;
                
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
            
            if(iteration % 50 == 0 && iteration > 0){
                for(int i = 0; i < nDestroy; i++)
                    destroyWeights[i] = 0.8 * destroyWeights[i] + 0.2 * (destroyScores[i] / destroyCounts[i]);
                for(int i = 0; i < nRepair; i++)
                    repairWeights[i] = 0.8 * repairWeights[i] + 0.2 * (repairScores[i] / repairCounts[i]);
            }
            
            T *= alpha;
            noImprove++;
            
            if(noImprove > 500){
                restartCount++;

                if(restartCount % 2 == 0){
                    cerr << "[Restart] Using Simplex solution.\n";
                    curAssignment = cachedSimplexSol;
                } else {
                    cerr << "[Restart] Using Greedy solution.\n";
                    curAssignment = greedySol;
                }

                curMerit = Prepare::meritFunction(curAssignment);
                noImprove = 0;
            }

            if(iteration % 1000 == 0){
                cerr << "Iteration " << iteration << ": Current Merit = " << curMerit << ", Best Merit = " << bestMerit << "\n";
            }
        }
        
        return bestAssignment;
    }
    
    static int rouletteSelect(const vector<double> &weights){
        double sum = 0;
        for(auto w : weights) sum += w;
        if(sum == 0) return rng() % weights.size();
        uniform_real_distribution<double> dist(0.0, sum);
        double r = dist(rng);
        double cum = 0;
        for(int i = 0; i < (int)weights.size(); i++){
            cum += weights[i];
            if(r <= cum) return i;
        }
        return weights.size() - 1;
    }
};
signed main(){

    // if(ifstream("input.txt")){
    //     freopen("input.txt","r",stdin);
    //     freopen("output.txt","w",stdout);
    // }

    Prepare::readData();
    Prepare::preProcess();
    
    cerr << "Running ALNS Algorithm...\n";
    vector<int> bestAssignment = ALNS::runAlgorithm(50000, 360.0);
    
    // Output results
    double bestMerit = Prepare::meritFunction(bestAssignment);
    cerr << Prepare::meritFunction(bestAssignment, true) << "\n";
    cerr << "Best Merit: " << bestMerit << "\n";
    cerr << "Solution: ";
    for(int i = 0; i < nNodes; i++){
        cerr << bestAssignment[i] << " ";
        cout << bestAssignment[i] << " ";
    }
    cerr << "\n";
    cout << "\n";

    
    return 0;
}