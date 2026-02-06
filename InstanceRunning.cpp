#include<bits/stdc++.h>

using namespace std;

void runAlgorithm(vector<pair<string,pair<int,int>>> testCases){
    for(int instance_id = 0;instance_id < (int)testCases.size();instance_id++){
        for(int id = testCases[instance_id].second.first;id <= testCases[instance_id].second.second;id++){
            string inputFile = "instances/" + testCases[instance_id].first + to_string(id) + ".dat";
            string outputFile = "results/" + testCases[instance_id].first + to_string(id) + ".out";
            system(("ProposedAlgorithm.exe < " + inputFile + " > " + outputFile).c_str());
            cout << "Processed " << inputFile << " -> " << outputFile << "\n";
        }

    }
}
void compareResults(vector<pair<string,pair<int,int>>> testCases){
    // Format current local date/time as YYYY-MM-DD HH:MM:SS
    time_t now = time(nullptr);
    tm* localNow = localtime(&now);
    char dateTimeBuf[20]{};
    strftime(dateTimeBuf, sizeof(dateTimeBuf), "%Y-%m-%d %H:%M:%S", localNow);
    string dateTimeString = dateTimeBuf;
    system(("echo " + dateTimeString + " > result.txt").c_str());

    for(int instance_id = 0;instance_id < (int)testCases.size();instance_id++){
        for(int id = testCases[instance_id].second.first;id <= testCases[instance_id].second.second;id++){
            string outputFile = "results\\" + testCases[instance_id].first + to_string(id) + ".out";
            string JuryFile = "AryResults\\" + testCases[instance_id].first + to_string(id) + ".out";
            string inputFile = "instances\\" + testCases[instance_id].first + to_string(id) + ".dat";
            string s = "#Instance " + testCases[instance_id].first + to_string(id) + ": ";
            system(("echo " + s + " >> result.txt").c_str());
            system(("copy " + inputFile + " checker.txt").c_str());
            system(("type " + JuryFile + " >> checker.txt").c_str());
            system(("type " + outputFile + " >> checker.txt").c_str());
            system("checker.exe < checker.txt >> result.txt");
            cout << "Compared results for " << outputFile << " and " << JuryFile << "\n";
            
        }
    }
}
signed main(){
    vector<pair<string,pair<int,int>>> testCases = {
        //{"Center486_G", {0, 9}},
        {"Center600_G", {0, 9}},
        {"Center726_G", {0, 9}},
        {"Corners486_G", {0, 9}},
        {"Corners600_G", {0, 9}},
        {"Corners726_G", {0, 9}},
        // {"Diagonal486_G", {0, 9}},
        // {"Diagonal600_G", {0, 9}},
        // {"Diagonal726_G", {0, 9}},
        // {"planar500_G", {0, 9}},
        // {"planar600_G", {0, 9}},
        // {"planar700_G", {0, 9}},
    };

    runAlgorithm(testCases);
    compareResults(testCases);

}