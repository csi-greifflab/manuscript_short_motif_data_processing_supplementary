#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <utility>
#include <string>


int A_coef = ('A' - 'C') * ('A' - 'G') * ('A' - 'T');
int C_coef = ('C' - 'A') * ('C' - 'G') * ('C' - 'T');
int G_coef = ('G' - 'A') * ('G' - 'C') * ('G' - 'T');
int T_coef = ('T' - 'A') * ('T' - 'C') * ('T' - 'G');
char reverse_complement(int x) {

    return (x - 'C') * (x - 'G') * (x - 'T') / A_coef * 'A' +
           (x - 'A') * (x - 'G') * (x - 'T') / C_coef * 'C' +
           (x - 'A') * (x - 'C') * (x - 'T') / G_coef * 'G' + 
           (x - 'A') * (x - 'C') * (x - 'G') / T_coef * 'T'; 
}

int main(int argc, char** argv) {
    if (argc != 3) { 
        std::cout << "usage:\n\tcount_unique_cdr3s.out input_filename output_filename" << std::endl;
        return 0;
    }

    std::map <std::string, std::string> codontable = {
        {"ATA","I"}, {"ATC","I"}, {"ATT","I"}, {"ATG","M"},
        {"ACA","T"}, {"ACC","T"}, {"ACG","T"}, {"ACT","T"},
        {"AAC","N"}, {"AAT","N"}, {"AAA","K"}, {"AAG","K"},
        {"AGC","S"}, {"AGT","S"}, {"AGA","R"}, {"AGG","R"},
        {"CTA","L"}, {"CTC","L"}, {"CTG","L"}, {"CTT","L"},
        {"CCA","P"}, {"CCC","P"}, {"CCG","P"}, {"CCT","P"},
        {"CAC","H"}, {"CAT","H"}, {"CAA","Q"}, {"CAG","Q"},
        {"CGA","R"}, {"CGC","R"}, {"CGG","R"}, {"CGT","R"},
        {"GTA","V"}, {"GTC","V"}, {"GTG","V"}, {"GTT","V"},
        {"GCA","A"}, {"GCC","A"}, {"GCG","A"}, {"GCT","A"},
        {"GAC","D"}, {"GAT","D"}, {"GAA","E"}, {"GAG","E"},
        {"GGA","G"}, {"GGC","G"}, {"GGG","G"}, {"GGT","G"},
        {"TCA","S"}, {"TCC","S"}, {"TCG","S"}, {"TCT","S"},
        {"TTC","F"}, {"TTT","F"}, {"TTA","L"}, {"TTG","L"},
        {"TAC","Y"}, {"TAT","Y"}, {"TAA","_"}, {"TAG","_"},
        {"TGC","C"}, {"TGT","C"}, {"TGA","_"}, {"TGG","W"},
    };

    std::map<std::int64_t, std::string> codon_id_to_aa;
    for (auto p : codontable) {
            std::string codon = p.first;
            std::string aa = p.second;

            std::int64_t c1 = codon[0] - 'A';
            std::int64_t c2 = codon[1] - 'A';
            std::int64_t c3 = codon[2] - 'A';
            std::int64_t codon_id = (c1 << 10) + (c2 << 5) + c3;
            codon_id_to_aa[codon_id] = aa;
        }

    std::vector<char> complement_nt(20);
    complement_nt[0] = 19;
    complement_nt[2] = 6;
    complement_nt[6] = 2;
    complement_nt[19] = 0;

    char* get_complement = &complement_nt[0];

    std::string in_filename = argv[1];
    std::string out_filename = argv[2];

    // std::ifstream infile("/storage/andreisl/data/for_brij/new_merged.fastq.assembled.fastq");
    // std::ofstream outfile("/storage/andreisl/data/for_brij/new_clean_sequences.fasta");
    std::ifstream infile(in_filename);
    std::ofstream outfile(out_filename);
    std::set<std::vector<std::int64_t>> nt_cdr3s;   
    std::map<std::vector<std::int64_t>, int> aa_cdr3s;   


    std::string line;
    bool cont = true;
    std::int64_t line_id = 0;
    std::int64_t wrong_length = 0;
    std::string header;
    int forward_found = 0;
    int reverse_found = 0;
    while (cont) {
        for (std::int64_t i = 0; i < 4 && cont; ++i) {
            cont = (bool) std::getline(infile, line);
            char* line_start = &line[0];
            if (!cont) {
                break;
            }

            if (i == 0) {
                header = line;
                header[0] = '>';
                continue;
            }
            if (i != 1) {
                continue;
            }
            if (line_id % 1000000 == 0) {
                std::cout << line_id / 1000000 << std::endl;
            }
            ++line_id;
            size_t len = line.size();
            if (len < 208) {
                ++wrong_length;
                continue;
            }

            bool rc = false;
            //GGATT -- rc_backward
            //ACCCG -- forward
            if ((*(line_start + 1) == 'G') && (*(line_start + 2) == 'G') && (*(line_start + 3) == 'A') && (*(line_start + 4) == 'T') && (*(line_start + 5) == 'T')) {
                rc = true;
            } else {
                if ((*(line_start + 1) != 'A') || (*(line_start + 2) != 'C') || (*(line_start + 3) != 'C') || (*(line_start + 4) != 'C') || (*(line_start + 5) != 'G')) {
                    continue;
                }
            }
            if (rc) {
                ++forward_found;
            } else {
                ++reverse_found;
            }

            std::vector<std::int64_t> aa_cdr3(52-39);
            std::vector<std::int64_t> nt_cdr3((52-39)*3);
            for (size_t ind = 39; ind < 52; ++ind) {
                std::int64_t c1;
                std::int64_t c2;
                std::int64_t c3;
                std::int64_t codon_id;
                if (!rc) {
                    c1 = *(line_start + ind * 3 + 1) - 'A';
                    c2 = *(line_start + ind * 3 + 2) - 'A';
                    c3 = *(line_start + ind * 3 + 3) - 'A';
                } else {
                    c1 = *(line_start + len - 1 - ind * 3 - 1) - 'A';
                    c1 = *(get_complement + c1);
                    c2 = *(line_start + len - 1 - ind * 3 - 2) - 'A';
                    c2 = *(get_complement + c2);
                    c3 = *(line_start + len - 1 - ind * 3 - 3) - 'A';
                    c3 = *(get_complement + c3);
                }
                if ((c1 != 0 && c1 != 2 && c1 != 6 && c1 != 19) || 
                    (c2 != 0 && c2 != 2 && c2 != 6 && c2 != 19) ||
                    (c3 != 0 && c3 != 2 && c3 != 6 && c3 != 19)) {
                    continue;
                }
                nt_cdr3[(ind - 39) * 3 + 0] = c1;
                nt_cdr3[(ind - 39) * 3 + 1] = c2;
                nt_cdr3[(ind - 39) * 3 + 2] = c3;
                codon_id = (c1 << 10) + (c2 << 5) + c3; 
                aa_cdr3[ind - 39] = codon_id;

            }

            nt_cdr3s.insert(nt_cdr3);
            if (aa_cdr3s.find(aa_cdr3) == aa_cdr3s.end()) {
                aa_cdr3s[aa_cdr3] = 0;
            }
            ++aa_cdr3s[aa_cdr3];
        }
    }

    outfile << "aa_cdr3_seq" << ',' << "count\n";
    for (auto const& p : aa_cdr3s) {
        std::string aa_cdr3_seq;
        for (int c : p.first) {
            aa_cdr3_seq += codon_id_to_aa[c];
        }
        outfile << aa_cdr3_seq << ',' << p.second << '\n';
    }
    outfile << std::endl;

    std::cout << line_id << std::endl;
    std::cout << "forward reads: " << forward_found << std::endl;
    std::cout << "rc reads: " << reverse_found << std::endl;
    std::cout << "total ok reads: " << forward_found + reverse_found << std::endl;
    std::cout << "unique nt cdr3s: " << nt_cdr3s.size() << std::endl;
    std::cout << "unique aa cdr3s: " << aa_cdr3s.size() << std::endl;
}
