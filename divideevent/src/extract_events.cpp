#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cstdint>
#include <sstream>

constexpr uint32_t MAGIC  = 0x45564E54; // 'EVNT'
constexpr size_t WORD_SIZE = 4;
constexpr size_t HEADER_WORDS = 8;

bool read_word(std::ifstream &f, uint32_t &w) {
    char bytes[WORD_SIZE];
    if (!f.read(bytes, WORD_SIZE)) return false;
    // Little Endian
    w = static_cast<uint32_t>(static_cast<unsigned char>(bytes[0]) |
                              (static_cast<unsigned char>(bytes[1]) << 8) |
                              (static_cast<unsigned char>(bytes[2]) << 16) |
                              (static_cast<unsigned char>(bytes[3]) << 24));
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 5) {
        std::cerr << "Usage: ./extract_events input.dat output_directory start_event end_event\n";
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_dir = argv[2];
    uint32_t start_event = std::stoul(argv[3]);
    uint32_t end_event   = std::stoul(argv[4]);

    std::ifstream fin(input_file, std::ios::binary);
    if (!fin) {
        std::cerr << "Failed to open input file: " << input_file << "\n";
        return 1;
    }

    // Output filename auto
    std::string filename = input_file.substr(input_file.find_last_of("/\\") + 1); // run00001copy.dat
    std::string base = filename.substr(0, filename.size() - 4); // run00001copy
    if (base.size() > 4 && base.substr(base.size()-4) == "copy") {
        base = base.substr(0, base.size()-4); // run00001
    }
    std::ostringstream filename_out;
    filename_out << output_dir << "/" << base << "_" << start_event << "_" << end_event << ".dat";
    std::ofstream fout(filename_out.str(), std::ios::binary);
    if (!fout) {
        std::cerr << "Failed to create output file\n";
        return 1;
    }
    std::cout << "Output file: " << filename_out.str() << std::endl;

    uint32_t word;
    std::vector<uint32_t> current_event;
    bool in_event = false;
    uint32_t current_event_num = 0;
    uint32_t max_event_num = 0;  // 最大イベント番号を記録
    std::vector<uint32_t> lookahead;

    while (read_word(fin, word)) {
        // Event header detection
        if (word == MAGIC) {
            // Look 8 words ahead
            std::streampos pos = fin.tellg();
            lookahead.clear();
            uint32_t tmp;
            bool valid = true;
            for (size_t k = 0; k < HEADER_WORDS; k++) {
                if (!read_word(fin, tmp)) {
                    valid = false;
                    break;
                }
                lookahead.push_back(tmp);
            }
            fin.seekg(pos); // rewind after lookahead
            if (valid && lookahead.size() >= HEADER_WORDS && lookahead[HEADER_WORDS - 1] == MAGIC) {
                // New event start
                if (in_event) {
                    // Write previous event if in range
                    if (current_event_num >= start_event && current_event_num <= end_event) {
                        fout.write(reinterpret_cast<char*>(current_event.data()), current_event.size() * WORD_SIZE);
                    }
                    current_event.clear();
                }
                in_event = true;

                // 2 words after MAGIC is the event number
                current_event_num = lookahead[1];
                if (current_event_num > max_event_num) {
                    max_event_num = current_event_num;  // update maximum
                }

                // Stop early if exceeded end_event
                if (current_event_num > end_event) break;
            }
        }

        // Accumulate event data
        if (in_event) {
            current_event.push_back(word);
        }
    }

    // Write last event if in range
    if (in_event && current_event_num >= start_event && current_event_num <= end_event) {
        fout.write(reinterpret_cast<char*>(current_event.data()), current_event.size() * WORD_SIZE);
    }

    fin.close();
    fout.close();

    // Check and warn if end_event exceeds actual max event
    std::cout << "Extraction complete.\n";
    std::cout << "Max event number in file: " << max_event_num << std::endl;
    if (end_event > max_event_num) {
        std::cerr << "Warning: Specified end_event (" << end_event
                  << ") exceeds the maximum event number in file (" << max_event_num << ").\n";
    }

    return 0;
}
