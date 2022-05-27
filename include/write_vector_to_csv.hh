
#include "includes_and_types.hh"

#pragma once

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

void
write_vector_to_csv
    ( const vector_t& vec
    , const std::string& file_path
    );