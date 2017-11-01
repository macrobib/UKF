#include "read_csv.h"
#include <string>
#include <fstream>
#include <sstream>

CSVHandler::CSVHandler(std::string file):csvfile(file), streamfile(file, std::fstream::in), 
                        init_status(false){
}

CSVHandler::~CSVHandler(){
    streamfile.close();
}


auto CSVHandler::split_string(const std::string& s, const char* delim) -> std::vector<std::string>{
    std::istringstream ss(s);
    std::string item;
    std::vector<std::string> output;
    while(std::getline(ss, item, *delim))
        output.push_back(std::move(item));
    return output;
}


auto CSVHandler::stream_handler(std::function<void(dstream_t&)> handler) -> bool{
    std::string line;
    auto status = true;
    static int count = 0;
    try{
        while(!streamfile.eof()){
            count++;
            std::getline(streamfile, line);
            auto output = split_string(line, "\t");
            std::string type = output.front();
            std::swap(output[0], output[output.size() - 1]);
            output.pop_back();
            data_capsule.type = type;
            data_capsule.px_rho = std::stod(output[1]);
            data_capsule.py_theta =  std::stod(output[2]);
            if(type == "R"){
                data_capsule.rhod =  std::stod(output[3]);
                data_capsule.timestamp =  std::stod(output[4]);
            }
            else{
                data_capsule.timestamp =  std::stod(output[3]);
                data_capsule.rhod = 0.0;
            }
            handler(data_capsule);
            if(count == 500)
                break;
        }
    }
    catch(std::ios_base::failure &fl){
        status = false;
        std::cout<< fl.what()<< std::endl;
    }
    return status;
}

auto CSVHandler::load_data_file(std::function<void(dstream_t&)> handler) -> void{
        //params.
        std::string type;
        double a, b, c, d, e, f;
        csvfile.read_header(io::ignore_extra_column,"type", "a", "b", "c", "d", "e", "f");
        while(csvfile.read_row(type, a, b, c, d, e, f)){
            data_capsule.type = type;
            data_capsule.px_rho = a;
            data_capsule.py_theta = b;
            if(type == "L"){
                data_capsule.timestamp = c;
                data_capsule.rhod = 0;
            }
            else{
                data_capsule.rhod = c;
                data_capsule.timestamp = d;
            }
            handler(data_capsule);
        }
}

//Get single data row at specific count.
auto CSVHandler::get_single_data_stream(int count) -> dstream_t{

    csvfile.read_header(io::ignore_extra_column,"type", "a", "b", "c", "d", "e", "f");
    std::string type;
    double a, b, c, d, e, f;
    if(csvfile.read_row(type, a, b, c, d, e, f)){
        data_capsule.type = type;
        data_capsule.px_rho = a;
        data_capsule.py_theta = b;
        if(type == "L"){
            data_capsule.timestamp = c;
            data_capsule.rhod = 0;
        }
        else{
            data_capsule.rhod = c;
            data_capsule.timestamp = d;
        }
    }
    return data_capsule;
}
