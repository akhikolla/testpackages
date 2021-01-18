
#include "readdf.h"

#include "memorymap.h"
#include "dataset.h"

#include <string>
#include <vector>
#include <map>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DataFrame readDF(
        String path,
        SEXP columnsReq,
        bool headerOnly)
{
    try
    {
        MemoryMap *mm = MemoryMap::attach(path);
        DataSet &dataset = *DataSet::retrieve(mm);

        int columnCount = dataset.columnCount();
        int rowCount = 0;
        int rowCountExFiltered = 0;

        if ( ! headerOnly)
        {
            rowCount = dataset.rowCount();
            rowCountExFiltered = dataset.rowCountExFiltered();
        }

        CharacterVector rowNames(rowCountExFiltered);

        int rowNo = 0;
        int colNo = 0;

        for (int i = 0; i < rowCount; i++)
        {
            if ( ! dataset.isRowFiltered(i))
                rowNames[rowNo++] = String(std::to_string(i+1));
        }

        bool readAllColumns;
        StringVector columnsRequired;
        List columns;
        CharacterVector columnNames;
        int offset = 0;

        if (Rf_isNull(columnsReq))
        {
            readAllColumns = true;

            // but not filters
            for (int i = 0; i < dataset.columnCount(); i++)
            {
                if (dataset[i].columnType() == ColumnType::FILTER)
                {
                    offset = i + 1;
                    columnCount--;
                }
                else
                {
                    break; // filters are only at the beginning of the dataset
                }
            }

            columns = List(columnCount);
            columnNames = CharacterVector(columnCount);
        }
        else
        {
            readAllColumns = false;
            columnsRequired = as<StringVector>(columnsReq);
            columns = List(columnsRequired.size());
            columnNames = CharacterVector(columnsRequired.size());
        }

        for (int i = 0; i < columnCount; i++)
        {
            Column column = dataset[offset + i];
            string columnName = column.name();

            if ( ! readAllColumns)
            {
                bool required = false;
                StringVector::iterator itr;
                for (itr = columnsRequired.begin(); itr != columnsRequired.end(); itr++)
                {
                    if (as<string>(*itr) == columnName)
                        required = true;
                }
                if ( ! required)
                    continue;
            }
            else
            {
                if (column.columnType() == ColumnType::FILTER)
                    continue;
            }

            columnNames[colNo] = String(columnName);

            if (column.columnType() == ColumnType::FILTER)
            {
                columns[colNo] = LogicalVector(rowCountExFiltered, true);
            }
            else if (column.dataType() == DataType::DECIMAL)
            {
                NumericVector v(rowCountExFiltered, NumericVector::get_na());
                rowNo = 0;

                for (int j = 0; j < rowCount; j++)
                {
                    if ( ! dataset.isRowFiltered(j))
                    {
                        if ( ! column.shouldTreatAsMissing(j))
                            v[rowNo] = column.raw<double>(j);
                        rowNo++;
                    }
                }

                columns[colNo] = v;
            }
            else if (column.dataType() == DataType::INTEGER && ! column.hasLevels())
            {
                IntegerVector v(rowCountExFiltered, IntegerVector::get_na());
                rowNo = 0;

                for (int j = 0; j < rowCount; j++)
                {
                    if ( ! dataset.isRowFiltered(j))
                    {
                        if ( ! column.shouldTreatAsMissing(j))
                            v[rowNo] = column.raw<int>(j);
                        rowNo++;
                    }
                }

                if (column.measureType() == MeasureType::ID)
                    v.attr("jmv-id") = true;

                columns[colNo] = v;
            }
            else if (column.dataType() == DataType::TEXT &&
                     column.measureType() == MeasureType::ID)
            {
                StringVector v(rowCountExFiltered, StringVector::get_na());
                rowNo = 0;

                for (int j = 0; j < rowCount; j++)
                {
                    if ( ! dataset.isRowFiltered(j))
                    {
                        if ( ! column.shouldTreatAsMissing(j))
                            v[rowNo] = String(column.raws(j));
                        rowNo++;
                    }
                }

                v.attr("jmv-id") = true;
                columns[colNo] = v;
            }
            else
            {
                int MISSING = IntegerVector::get_na();

                // populate levels

                vector<LevelData> m = column.levels();

                int nLevels = column.levelCountExFilteredExMissing();
                CharacterVector levels = CharacterVector(nLevels);
                IntegerVector values = IntegerVector(nLevels);

                map<int, int> indexes;
                int jli = 0;
                int rli = 0;

                vector<LevelData>::iterator itr = m.begin();
                for (; itr != m.end(); itr++)
                {
                    LevelData &p = *itr;
                    if (p.filtered() == false && p.treatAsMissing() == false)
                    {
                        int value;

                        if (column.dataType() == DataType::TEXT)
                            value = jli;
                        else
                            value = p.ivalue();

                        indexes[value] = rli + 1;
                        values[rli] = value;
                        levels[rli] = String(p.label());

                        jli++;
                        rli++;
                    }
                    else
                    {
                        jli++;
                    }

                }

                // populate cells

                IntegerVector v(rowCountExFiltered, MISSING);
                rowNo = 0;

                for (int j = 0; j < rowCount; j++)
                {
                    if ( ! dataset.isRowFiltered(j))
                    {
                        int value = column.raw<int>(j);
                        if (value != INT_MIN)
                        {
                            if ( ! column.shouldTreatAsMissing(j))
                                v[rowNo] = indexes[value];
                            else
                                v[rowNo] = MISSING;
                        }
                        rowNo++;
                    }
                }

                // assign levels

                v.attr("levels") = levels;

                if (column.dataType() == DataType::TEXT)
                {
                    if (column.measureType() == MeasureType::ORDINAL)
                        v.attr("class") = CharacterVector::create("ordered", "factor");
                    else
                        v.attr("class") = "factor";
                }
                else
                {
                    v.attr("values") = values;

                    if (column.measureType() == MeasureType::ORDINAL)
                        v.attr("class") = CharacterVector::create("ordered", "factor");
                    else
                        v.attr("class") = CharacterVector::create("factor");
                }

                if ( ! column.trimLevels() && column.hasUnusedLevels())
                    v.attr("jmv-unused-levels") = true;

                if (column.measureType() == MeasureType::ID)
                    v.attr("jmv-id") = true;

                columns[colNo] = v;
            }

            colNo++;
        }

        if (colNo < columnsRequired.size())
        {
            columns.erase(colNo, columnsRequired.size());
            columnNames.erase(colNo, columnsRequired.size());
        }

        columns.attr("names") = columnNames;
        columns.attr("row.names") = rowNames;
        columns.attr("class") = "data.frame";

        return columns;
    }
    catch (std::exception &ex) {
        forward_exception_to_r(ex);
    }
    catch (...) {
        ::Rf_error("c++ exception (unknown reason)");
    }

    return DataFrame();
}
