
#include "write_vector_to_csv.hh"

void
write_vector_to_csv
    ( const vector_t& vec
    , const std::string& file_path
    )
{
    std::ofstream file(file_path.c_str());
    file << vec.format(CSVFormat);
} 