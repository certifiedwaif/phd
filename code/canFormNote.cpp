// canFormNote.cpp
#include <iostream>
#include <string>

using namespace std;

bool canFormNote(string poster, string note)
{
	int letter_counts[26] = {};
	
	for (char c : poster) {
		letter_counts[c - 'A']++;
	}

	for (char c : note) {
		letter_counts[c - 'A']--;
		if (letter_counts[c - 'A'] < 0)
			return false;
	}

	return true;
}

int main(int argc, char **argv)
{
	cout << (canFormNote("AB", "A") ? "true" : "false") << endl;

	return 0;
}