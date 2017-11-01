#include "csv.h"
#include <iostream>
#include <fstream>
#include <string>

typedef struct{
        std::string type; //Data type
        double px_rho;   //X position.
        double py_theta;   //Y Position.
        double rho;  //Bearing.
        double rhod; //Bearing rate.
        double timestamp; //Current time stamp.
    }dstream_t;

class CSVHandler{
public:
    CSVHandler();
    ~CSVHandler();
    CSVHandler(std::string file);
    auto load_data_file(std::function<void(dstream_t&)> handler) -> void;
    auto get_single_data_stream(int count = 0) -> dstream_t;
    auto stream_handler(std::function<void(dstream_t&)>) -> bool;
private:
    io::CSVReader<7> csvfile;
    std::fstream streamfile;
    bool init_status;
    dstream_t data_capsule;

    auto split_string(const std::string&, const char*) -> std::vector<std::string>;
};
