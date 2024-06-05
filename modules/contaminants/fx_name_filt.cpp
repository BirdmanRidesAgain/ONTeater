#include <iostream>
#include <fstream>
#include <vector>

/* 
Usage: ./fx_name_filt test_longread_mtDNA.fq output.csv | gzip > test_longread_mtDNA_filt.fq.gz
    OR
Usage: zcat test_longread_mtDNA.fq.gz | ./fx_name_filt stdin output.csv | gzip > test_longread_mtDNA_filt.fq.gz

 Compile command is: g++ fx_name_filt.cpp -o fx_name_filt
    OR
*/

/*
NEEDED FIXES:
    get it to accept fasta files
*/

// HELPER FUNCTIONS
std::string get_filetype(std::string fastx, std::string delimiter = "." ) {
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
std::vector<std::string> parse_queries(std::string query) { //return a vector of taxids
    std::ifstream input(query);
    std::vector<std::string> queries;
    std::string line;

    //check for okay file
    if (!input.good()) {
        std::cerr << "Error opening '" << query << "'. Exiting program." << std::endl;
        queries.push_back("-99");
        return queries;
    }

    while(std::getline(input, line).good()) {
        //do stuff while the line isn't empty
        queries.push_back(line);
    }
    return(queries);
}
bool filter_on_name(std::string line, std::vector<std::string> queries) {
    bool isFound = false;
    for (std::vector<std::string>::iterator i = queries.begin(); i != queries.end(); i++) {
        if (line.find(*i) != std::string::npos && line.find(*i) + (*i).length() == line.length()) {
            isFound = true;
        }
    }
    return(!isFound);
}

void parse_fasta(std::string fastx, std::vector<std::string> queries) {
    std::cout << "Not yet implemented. Try later." << std::endl;
}
void parse_fastq(std::string fastx, std::vector<std::string> queries) {
    std::ifstream input(fastx);
    if (!input.good()) {
        std::cerr << "Error opening '" << fastx << "'. Exiting program." << std::endl;
    }
    std::string line, name, content, qual; 
    int counter = 0; 

    // begin reading fastq
    while (std::getline(input, line).good()) {
        if (!line.empty()) {
            if ((counter % 4) == 0 ) {
                if ( !name.empty() && !content.empty() && !qual.empty()) {
                    bool queryPass = filter_on_name(name, queries);
                    if (queryPass) {
                        std::cout << name << std::endl << content << std::endl << "+" << std::endl << qual << std::endl;
                    }
                }
                name = line;
            }
            else if((counter % 4 == 1)) {
                content = line;
            }
            else if ((counter % 4 == 3)) {
                qual = line;
            }
            counter++;
        }
    }
}



// DRIVER CODE
int main( int argc, char **argv ) {
    // 1. PARSE ARGUMENTS
        //If no arguments, output usage:
        if (argc == 0 ) {
        std::cerr << "Usage: " << argv[0] << " [infile.fx] [kraken2_report.txt]" << std::endl;
        std::cerr << "  To accept data from stdin, set 'infile.fx' to 'stdin'. You still must supply an arg for kraken2_report.txt." << std::endl;
        return -1;
    }

    std::string fastx = argv[1]; // the fastx file, if not provided in stdin.
    std::string query = argv[2]; // the query file. MUST be provided; if absent no filtering

    // 2. READ IN QUERY FILE
    std::vector<std::string> queries = parse_queries(query);
    if (queries[0] == "-99") { return -1; } //invalid queries


    // 3. TAKE IN FASTX FROM STDIN IF APPLICABLE
        //BROKEN NEEDS FIXING
        //MAYBE PUT IN OWN FUNCTION?
    if (fastx == "stdin") {
        //Program will hang if no stdin is provided.
            //check for stdin output somehow?
        //Program only accepts fq input from stdin - fa will cause malfunctions.
        //STDIN DOES NOT WORK
        std::string line, name, content, qual; 
        int counter = 0; 

        // begin reading stdin
        while (std::getline(std::cin, line).good()) {
            if (!line.empty()) {
                if ((counter % 4) == 0) {
                    if (!name.empty() && !content.empty() && !qual.empty()) {
                        bool queryPass = filter_on_name(name, queries);
                        if (queryPass) {
                            std::cout << name << std::endl << content << std::endl << "+" << std::endl << qual << std::endl;
                        }
                    }
                    name = line;
                }
                else if ((counter % 4) == 1) {
                    content = line;
                }
                else if ((counter % 4) == 3) {
                    qual = line;
                }
                counter++;
            }
        }
    }

    // 4. TAKE IN FASTX FROM FILE IF APPLICABLE
    else {
        std::string filetype = get_filetype(fastx);
    // Read in fastx
    if ( filetype == "fa" ) {
        //treat input like fasta
        parse_fasta(fastx, queries);
    }
    else if ( filetype == "fq") {
        //treat input like fastq
        parse_fastq(fastx, queries);
    } else {
        std::cout << "Filetype unknown. Exiting program." << std::endl;
        return -1;
    }
    }
    return 0;
}
