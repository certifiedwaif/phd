// canFormNote.cpp
#include <iostream>
#include <string>
#include <vector>
#include <utility>

using namespace std;

bool canFormNote(string& poster, string& note)
{
	int letter_counts[26] = {};
	
	for (char c : poster) {
		letter_counts[c - 'A'] += 1;
	}

	for (char c : note) {
		letter_counts[c - 'A'] -= 1;
		if (letter_counts[c - 'A'] < 0)
			return false;
	}

	return true;
}

string showBool(bool b)
{
	return b ? "true" : "false";
}

int main()
{
	vector< pair< pair<string, string>, bool > > v = {
		{{"AB", "A"}, true},
		{{"A", "AB"}, false},
		{{"", "A"}, false},
		{{"A", ""}, true},
	};

	for (auto& p : v) {
		auto& arg_pair = p.first;
		string& poster = arg_pair.first;
		string& note = arg_pair.second;
		bool actual = canFormNote(poster, note);
		bool expected = p.second;

		cout << showBool(actual) << " " << showBool(expected) << endl;
	}

	return 0;
}