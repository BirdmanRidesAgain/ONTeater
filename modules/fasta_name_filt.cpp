#include <iostream>
#include <fstream>
using std::string;

//program designed to read in a fasta, classify it based on the header line, and print it back out
/* IMPROVEMENTS:
Allow it to take gzipped fq files
allow it to take files from command line
allow it to filter on including everything as well as removing all queries
allow it to filter based on string content
let it use fasta files
let it write to a designated output file

rework the loop so it stops reading the next three lines if the name is filtered out
*/
//kraken:taxid|562

bool filter_on_substring(string line, string query = "") {
    //take query argument from input args and then dump all string names that include it
    //set invert flag to true if you want to filter out all names that DON'T have a particular string
    bool isFound = line.find(query) != string::npos;
    return (!isFound); // if 'found' is true, you fail the filter
}

int parse_fastq(string fastx, string query = "") {
    std::ifstream input(fastx); //create input file stream called "input" and associate it with arg. 1
    if (!input.good()) { // if the input is invalid, print error and return
        std::cerr << "Error opening '" << fastx << "'. Exiting program." << std::endl;
        return -1;
    }
    //initialize four strings called line, name and content
    string line; //'working' portion of fx file. either written to 'name' or content
    string name; //name - overwritten every time a new '>' line is read in
    string content; //all lines that aren't empty, begin with '>', or have ' '
    string qual; // phred-scaled quality score

    //Read in lines. This implementation requires that fastq files be written in single-line format
    int counter = 0; //line counter for fq file
    while( std::getline( input, line ).good() ) {
        //if the line isn't empty, do stuff
        if ( !line.empty() ) { //This check means that accidental blank lines don't screw up your parse
            if ( (counter % 4) == 0 ) { // line is a name line
                if ( !name.empty() && !content.empty() && !qual.empty()) {
                //FIXME if all these are filled, make filter decision 
                //FIXME change filter location so we don't waste time reading in lines that we know don't fit based on the title
                    bool filterPass = filter_on_substring(name, query); 
                    if (filterPass) {
                        std::cout << name << std::endl;
                        std::cout << content << std::endl;
                        std::cout << "+" << std::endl;
                        std::cout << qual << std::endl;
                    }
                }
                name = line;
            } else if ( (counter % 4) == 1 ) { //line is a content line
                content = line;
            } 
            // if counter is 2, it's a separator and we don't care about it
            else if ( ( counter % 4 == 3) ) { // line is a qual string
                qual = line;
            }
            counter++; //only activates if you've parsed an actual fastq line
        }
    }

    return 0;
}
// FIXME - parse_fasta does not filter anything or make use of query
int parse_fasta(string fastx, string query) {
        std::ifstream input(fastx); //create input file stream called "input" and associate it with arg. 1
    if (!input.good()) { // if the input is invalid, print error and return
        std::cerr << "Error opening '" << fastx << "'. Exiting program." << std::endl;
        return -1;
    }

    //initialize three strings called line, name and content
    std::string line; //'working' portion of fx file. either written to 'name' or content
    std::string name; //name - overwritten every time a new '>' line is read in
    std::string content; //all lines that aren't empty, begin with '>', or have ' '
    
    while( std::getline( input, line ).good() ) { //loop repeats for every line in file
        if ( line.empty() || line[0] == '>' ) { //initiates new fasta when 'line' is empty or a new name
            if ( !name.empty() ) { //if you just initiated a new fasta AND name is defined, print the fa
                std::cout << name << " : " << content << std::endl;
                name.clear(); //clear name from buffer
            }
            if ( !line.empty() ) { //will only activate when a name (beginning with '>')
                name = line.substr(1); //strips the '>' from the name
            }
            content.clear(); //clears content, so you don't concatenate fa's together
        } 
        else if ( !name.empty() ) {
            //if whitespace is found before the final space in the string, the fa is invalid
                // this means you need to clear name and content, to reset the cycle
            if ( line.find(' ') != std::string::npos ) {
                name.clear();
                content.clear();
            } else { //if line IS valid, cat to content
                content += line;
            }
        }
    }
    if ( !name.empty() ) { // Print out what was read from the last entry
        std::cout << name << " : " << content << std::endl;
    }
    return 0;
}
string get_filetype( string fastx, string delimiter = "." ) {
    std::string filetype;
    std::string ext = fastx.substr(fastx.find(delimiter) + 1);
    if ( ext == "fasta" || ext == "fa" || ext == "fasta.gz" || ext == "fa.gz" ) {
        filetype = "fa";
    }
    if ( ext == "fastq" || ext == "fq" || ext == "fastq.gz" || ext == "fq.gz") {
        filetype = "fq";
    } else { filetype = "unknown"; } 
    return filetype;
}

int main( int argc, char **argv ) {
    if (argc <= 1 ) { //if there are no arguments supplied to the invocation, print usage to error and return
        std::cerr << "Usage: " << argv[0] << " [infile]" << std::endl;
        return -1;
    }
    std::ifstream input(argv[1]); //create input file stream called "input" and associate it with arg. 1
    if (!input.good()) { // if the input is invalid, print error and return
        std::cerr << "Error opening '" << argv[1] << "'. Exiting program." << std::endl;
        return -1;
    }

    //if the file is a fastq, activate this one
    //std::string fasta = "test_mtDNA_contigs.fasta";
    string fastx = argv[1];
    string query = argv[2];

    //attempt to parse what kind of fastx file this is
    string filetype = get_filetype(fastx);
    if ( filetype == "fa" ) {
        int fasta_success = parse_fasta(fastx, query);
    }
    else if ( filetype == "fq") {
        int fastq_success = parse_fastq(fastx, query);
    } else {
        std::cout << "Filetype unknown. Exiting program." << std::endl;
        return -1;
    }
    return 0;
}