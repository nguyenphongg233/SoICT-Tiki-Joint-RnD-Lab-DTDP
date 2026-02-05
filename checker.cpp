#include<bits/stdc++.h>

using namespace std;

#define read() ios_base::sync_with_stdio(0);cin.tie(0);cout.tie(0)
#define day() time_t now = time(0);char* x = ctime(&now);cerr<<"Right Now Is : "<<x<<"\n"

#define int long long
#define ii pair<int,int>
#define X first1
#define Y second 

const long long MAX = (int)1000 + 5;
const long long INF = (int)1e9;
const long long MOD = (int)1e9 + 7;
const long long N = 105;

int n,m,p;
struct Node{
    int id;
    double x,y;
    double w[3];
}node[MAX];
vector<pair<int,int>> adj[MAX];
double tau[] = {0.05,0.05,0.05};
vector<int> jury_assignment;
vector<int> user_assignment;

int dist[MAX][MAX];
int max_dist = 0;
void dijkstra(){
    for(int i = 0;i < n;i++){
        for(int j = 0;j < n;j++){
            dist[i][j] = INF;
        }
    }
    for(int s = 0;s < n;s++){
        priority_queue<ii,vector<ii>,greater<ii>> pq;
        dist[s][s] = 0;
        pq.push({0,s});
        while(!pq.empty()){
            ii top = pq.top();
            pq.pop();
            int d = top.first;
            int u = top.second;
            if(d > dist[s][u]) continue;
            for(auto edge : adj[u]){
                int v = edge.first;
                int w = edge.second;
                if(dist[s][v] > dist[s][u] + w){
                    dist[s][v] = dist[s][u] + w;
                    pq.push({dist[s][v],v});
                }
            }
        }
    }
    for(int i = 0;i < n;i++){
        for(int j = 0;j < n;j++){
            if(dist[i][j] < INF && dist[i][j] > max_dist){
                max_dist = dist[i][j];
            }
        }
    }
}

pair<long double, long double> calc(vector<int> &assignment){
    long double compactness_cost = 0.0;
    long double balance_cost = 0.0;
    vector<vector<int>> district(p + 5,vector<int>());
    long double w_sum[3] = {0,0,0};

    for(int i = 0;i < n;i++){
        district[assignment[i]].push_back(i);
        w_sum[0] += node[i].w[0];
        w_sum[1] += node[i].w[1];
        w_sum[2] += node[i].w[2];
    }

    w_sum[0] /= 1.0 * p;
    w_sum[1] /= 1.0 * p;
    w_sum[2] /= 1.0 * p;

    for(int i = 0;i < p;i++){
        // Balance Cost
        long double curr_w[3] = {0,0,0};
        for(auto u : district[i]){
            curr_w[0] += node[u].w[0];
            curr_w[1] += node[u].w[1];
            curr_w[2] += node[u].w[2];
            for(auto v : district[i]){
                compactness_cost = max(compactness_cost, (long double)dist[u][v] / max_dist);
            }
        }
        balance_cost += max({0.0L,(1 - tau[0]) * w_sum[0] - curr_w[0],curr_w[0] - (1 + tau[0]) * w_sum[0]}) / (double)w_sum[0];
        balance_cost += max({0.0L,(1 - tau[1]) * w_sum[1] - curr_w[1],curr_w[1] - (1 + tau[1]) * w_sum[1]}) / (double)w_sum[1];
        balance_cost += max({0.0L,(1 - tau[2]) * w_sum[2] - curr_w[2],curr_w[2] - (1 + tau[2]) * w_sum[2]}) / (double)w_sum[2];
    }
    return {compactness_cost, balance_cost};
}
signed main(){
	
	read();

    freopen("checker.txt","r",stdin);

    cin >> n;
    for(int i = 0;i < n;i++){
        cin >> node[i].id >> node[i].x >> node[i].y;
        for(int j = 0;j < 3;j++){
            cin >> node[i].w[j];
        }
    }
    cin >> m;
    for(int i = 0;i < m;i++){
        int u,v,w;
        cin >> u >> v >> w;
        adj[u].push_back({v,w});
        adj[v].push_back({u,w});
    }
    cin >> p >> tau[0] >> tau[1] >> tau[2];
    dijkstra();

    // Jury solution
    for(int i = 0,c;i < n;i++){
        if(cin >> c){
            jury_assignment.push_back(c);
        }  else{
            printf("%d Invalid solution from jury\n", 2);
            return 0;
        }
    }
    // User solution
    for(int i = 0,c;i < n;i++){
        if(cin >> c){
            user_assignment.push_back(c);
        }else{
            printf("%d Invalid solution from participant\n", -2);
            return 0;
        }
    }

    pair<long double, long double> jury_cost = calc(jury_assignment);
    pair<long double, long double> user_cost = calc(user_assignment);

    cout << fixed << setprecision(9) << "Jury Compactness: " << jury_cost.first << ", Balance: " << jury_cost.second << "\n";
    cout << fixed << setprecision(9) << "Participant Compactness: " << user_cost.first << ", Balance: " << user_cost.second << "\n";
    
    if(jury_cost.second < 1e-9){
        if(user_cost.second < 1e-9){
            // both valid
            if(user_cost.first < jury_cost.first - 1e-9){
                printf("%d Better solution! Jury_compactness = %Lf Participant_compactness = %Lf\n", 1, jury_cost.first, user_cost.first);
                return 0;
            }else if(fabs(user_cost.first - jury_cost.first) <= 1e-9){
                printf("%d Equal solution! Jury_compactness = %Lf Participant_compactness = %Lf\n", 0, jury_cost.first, user_cost.first);
                return 0;
            }else {
                printf("%d Worse solution! Jury_compactness = %Lf Participant_compactness = %Lf\n", -1, jury_cost.first, user_cost.first);
                return 0;
            }
        }else{
            // jury valid, user invalid
            printf("%d Balance violation solution from participant! Participant_balance = %Lf\n", -1, user_cost.second);
            return 0;
        }
    }else{
        if(user_cost.second < 1e-9){
            // jury invalid, user valid
            printf("%d Balance violation solution for Jury, Valid solution from participant! Participant_compactness = %Lf\n", 2, user_cost.first);
            return 0;
        }else{
            // both invalid
            if(user_cost.second + 1e-9 < jury_cost.second){
                printf("%d Less violation! Jury_balance = %Lf Participant_balance = %Lf\n", 1, jury_cost.second, user_cost.second);
                return 0;
            }else if(fabs(user_cost.second - jury_cost.second) <= 1e-9){
                printf("%d Equal violation! Jury_balance = %Lf Participant_balance = %Lf\n", 0, jury_cost.second, user_cost.second);
                return 0;
            }else {
                printf("%d More violation! Jury_balance = %Lf Participant_balance = %Lf\n", -1, jury_cost.second, user_cost.second);
                return 0;
            }
        }
    }
}

