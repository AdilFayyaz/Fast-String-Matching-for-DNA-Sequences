/*
	Compiled and Run on Visual Studio 2019
	M. Adil Fayyaz
	Dated: 15th May 2020
*/

#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include <iterator>
#include <map>
using namespace std;
typedef vector<vector<int> > TwoDMatrix;
typedef vector<int> OneDMatrix;
typedef vector<float> OneDFMatrix;
string inputfile;
string outputfile;

void PreProc_Calc(string P, string T, string E) {
	//Opening the Output File
	fstream out;
	out.open(outputfile.c_str(),ios::out);
	if (!out) {
		cout << "Output file didn't open!" << endl;
		exit(0);
	}
	//Initializations
	int m, range, n;
	m = P.length(); n = T.length(); range = E.length();
	TwoDMatrix shift(((m + 1) * range), vector<int>(((m + 1) * range), 1));
	OneDMatrix scan; OneDMatrix U; OneDMatrix safe;
	OneDFMatrix avr_shift; OneDFMatrix Ft(range, 0);
	scan.resize(m + 1); safe.resize(m + 1); avr_shift.resize(m + 1);
	Ft.resize(range);
	map <char, int> Sigma;
	for (int i = 0; i < range; i++) {
		Sigma.insert(pair<char, int>(E[i], i));
	}
	for (int i = 0; i < m; i++) {
		U.push_back(i + 1);
		safe[i] = 0;
	}
	cout << "Value of m is: " << m << endl;
	cout << "Probabilities: " << endl;
	//Calculate the probabilities
	for (int i = 0; i < range; ++i) {
		for (int x = 0; x < n; x++) {
			if (T[x] == E[i]) {
				Ft[i]++;
			}
		}
		Ft[i] /= n;
		cout << E[i]<<" "<<Ft[i] << endl;
	}
	cout << endl;

	//Preprocessing Stage
	for (int i = 1; i <= m; ++i) {
		vector<int>::iterator l2;
		for (l2 = U.begin(); l2 != U.end(); ++l2) {
			int l = *l2;
			avr_shift[l] = 0;
			for (map<char, int>::iterator s = Sigma.begin(); s != Sigma.end(); ++s) {
				for (int k = shift[l][s->second]; k <= m; ++k) {
					if (safe[k] == 0 and ((l - k) >= 0) and s->first == P[l - k]) {
						shift[l][s->second] = k; break;
					}
				}
				avr_shift[l] += Ft[s->second] * shift[l][s->second];
			}
		}
		//Finding the largest value l in U such that avr_shift[l] is maximum
		float max = INT_MIN;
		int val = 0;
		int iter = 0;
		vector<int>::iterator l;
		for (l = U.begin(); l != U.end(); ++l) {
			if (avr_shift[*l] > max) {
				max = avr_shift[*l];
				iter = *l;
			}
		}
		vector<int>::iterator J = find(U.begin(), U.end(), iter);
		if (J != U.end()) {
			int j = *J;
			U.erase(J);
			scan[i] = j - 1;
			for (int k = 1; k <= j - 1; ++k) {
				if (P[j - 1] != P[j - k - 1]) {
					safe[k - 1] = 1;
				}
			}
		}

	}

	//String Matching Algorithm 
	int w = 0;
	while (w <= n - m ) {
		int i = 1;
		while (i <= m && T[w + scan[i]] == P[scan[i]]) {
			++i;
		}
		if (i > m ) {
			cout << "Pattern found at index " << w << endl;
			if (out) {
				out << "Pattern found at index " << w << endl;
			}
			map<char, int>::iterator itr = Sigma.find(T[w + scan[m]]);
			w += shift[scan[m]][itr->second];
		}
		else {
			map<char, int>::iterator itr = Sigma.find(T[w + scan[i]]);
			w += shift[scan[i]][itr->second];
		}
	}
	out << endl << endl;
	out << "The Index position stated above gives the starting index from where the pattern '"<<P<<"' was found in the Text string";
	out.close();
}
int main() {
	string E = "ACGT";
	/*
		The pattern E (sigma) is the same for all test cases, because
		DNA sequences can be composed of combinations of only these 4 alphabets.
	*/
	while (1) {
		cout << "Enter the Input file's name: ";
		cin >> inputfile;
		ifstream in;
		in.open(inputfile.c_str(), ios::in);
		if (!in) {
			cout << "Input file does not exist\n";
			exit(0);
		}
		cout << endl;
		int ofn;
		cout << "Enter the Output file number: ";
		cin >> ofn;
		cout << endl;
		outputfile = "Output";
		outputfile += to_string(ofn);
		outputfile += ".txt";
		cout << "Output file is: " << outputfile << endl;
		
		string P; string T;
		if (!in) {
			cout << "Error Opening the file\n";
			exit(1);
		}
		in >> P;
		in >> T;
		PreProc_Calc(P, T, E);
		in.close();	
		char answer;
		cout << "Do you wish to read another file? (Y/N) ";
		cin >> answer;
		if (answer == 'N' || answer == 'n') { break; }
	}
	return 0;
}
/*
	Run Time Analysis: 
	Time for Preprocessing in BigO notation: O(m^2 *|sigma|)
	Time for String Matching Algorithm: O(m*n)
*/