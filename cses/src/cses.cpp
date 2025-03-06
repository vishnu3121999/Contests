#include <bits/stdc++.h>
using namespace std;

const int MOD = 1000000007;
unordered_map<long long, unordered_map<long long, long long>> memo;

long long power(long long a, long long b) {
    // Check if result is already computed
    if (memo.count(a) && memo[a].count(b)) {
        return memo[a][b];
    }

    long long result = 1;
    long long base = a % MOD;

    while (b > 0) {
        // If the current bit of b is set, multiply the result by base
        if (b & 1) {
            result = (result * base) % MOD;
        }
        // Square the base and reduce modulo MOD
        base = (base * base) % MOD;
        b >>= 1;
    }

    // Store result in memoization map
    memo[a][b] = result;
    return result;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(nullptr);

    int T;
    cin >> T;
    while (T-- > 0) {
        int a, b;
        cin >> a >> b;
        cout << power(a, b) << "\n";
    }

    return 0;
}
