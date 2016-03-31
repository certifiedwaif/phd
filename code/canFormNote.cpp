// canFormNote.cpp
#include <iostream>
#include <string>
#include <vector>
#include <utility>

using namespace std;

const bool canFormNote(const string& poster, const string& note)
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
	const vector< pair< pair<string, string>, bool > > v = {
		{{"AB", "A"}, true},
		{{"A", "AB"}, false},
		{{"", "A"}, false},
		{{"A", ""}, true},
	};

	for (const auto& p : v) {
		const auto& arg_pair = p.first;
		const string& poster = arg_pair.first;
		const string& note = arg_pair.second;
		const bool actual = canFormNote(poster, note);
		const bool expected = p.second;

		cout << showBool(actual) << " " << showBool(expected) << endl;
	}

	return 0;
}