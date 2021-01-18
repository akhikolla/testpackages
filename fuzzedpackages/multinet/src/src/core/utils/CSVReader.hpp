/**
 * History:
 * - 2018.01.01 file imported from version 1.0 of the multinet library
 */

#ifndef UU_CORE_UTILS_CSV_H_
#define UU_CORE_UTILS_CSV_H_

//#include <unordered_set>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <set>

namespace uu {
namespace core {

/**
 * Utility class to read a CSV file.
 */
class CSVReader
{
  private:
    std::ifstream infile_;
    std::string next_;
    std::string current_;
    bool has_next_;
    bool remove_trailing_spaces_;
    int next_row_number_;
    int current_row_number_;
    char field_separator_;
    char quote_;
    std::string quote_as_string_;
    std::string escaped_quote_;
    std::string comment_start_;
    size_t expected_num_fields_ = 0;
    size_t lines_skipped_ = 0;

  public:
    CSVReader();
    virtual
    ~CSVReader();

    /**
     * Opens a csv file, and sets the CSVReader before the first line.
     * @param path csv file
     */
    void
    open(
        const std::string& path
    );

    /**
     * Closes a csv file.
     */
    void
    close(
    );

    /**
     * Checks if the end of file has been reached.
     * @return true if there is at least one row still to be read
     */
    bool
    has_next(
    ) const;

    /**
     * Returns the next row from the file.
     * The CSVReader advances to the next line in the file.
     * Empty lines are automatically skipped.
     * Notice that this line can be pre-fetched during the previous reading operation. Therefore,
     * applying new settings (e.g., setting a new format for comments) may not apply immediately
     * and require one more call to get_next or get_next_raw_line.
     * @return a vector of strings, each representing a field
     */
    std::vector<std::string>
    get_next(
    );

    /**
     * Returns the next row from the file, without splitting it into fields.
     * The CSVReader advances to the next line in the file.
     * Empty lines are automatically skipped.
     * Notice that this line can be pre-fetched during the previous reading operation. Therefore,
     * applying new settings (e.g., setting a new format for comments) may not apply immediately
     * and require one more call to get_next or get_next_raw_line.
     * @return a string corresponding to the next line
     */
    std::string
    get_next_raw_line(
    );


    /**
     * Returns the current row from the file, without splitting it into fields.
     * The CSVReader remains at the same position (line) in the file.
     * Empty lines are automatically skipped.
     * @return a string corresponding to the current line
     * @throw OperationNotPermittedException if no get_next has been called yet
     */
    std::string
    get_current_raw_line(
    ) const;

    /**
     * @return the current row number.
     */
    int
    row_num(
    ) const;

    /**
     * Controls whether trailing blank characters at the end of a field should
     * be kept or removed.
     * The default value is false.
     * @param trim if true, the strings are trimmed
     */
    void
    trim_fields(
        bool trim
    );

    /**
     * Sets the character used to separate fields.
     * The default value is ','.
     * @param separator a character separating strings
     */
    void
    set_field_separator(
        char separator
    );

    /**
     * Sets the expected number of fields.
     * Lines with a different number of fields will be skipped.
     * The number of such lines can be retrieved using the skipped_lines() function.
     */
    void
    set_expected_num_fields(
        size_t expected_num_fields
    );


    /**
     * Returns the number of lines skipped because of an unexpected number of fields.
     */
    size_t
    skipped_lines(
    ) const;

    /**
     * Sets the characters that indicate a comment when at the beginning of
     * a line.
     * By default, comments are not handled.
     * If an empty string is passed to this method, comments are not handled.
     * @param comment_start indicates that the line should not be processed.
     */
    void
    set_comment(
        const std::string& comment_start
    );

    /**
     * Sets the character that indicates the start and end of a quoted field.
     * A quoted field starts when a field separator or the beginning of a line
     * is followed by the quote character.
     * Inside a quoted field, all consecutive occurrences of two quote characters
     * will be read as a single quote character, and occurrences of a field separator
     * will be considered as normal text and not break the field.
     * A quoted field ends when a new quote character, which is not the second in a sequence
     * of two as described above, is followed by a field separator character or the end of line.
     * Notice that single quote charaters not followed by a field separator character are
     * also allowed in a quoted field, and are considered as normal characters.
     *
     * Examples:
     * - "a, b, c" (followed by a comma or an end of line) is a single field.
     * - "a, b, ""c""" (followed by a comma or an end of line) is a single field, read as a, b, "c".
     * - "a, b, " c" (followed by a comma or an end of line) is also a single field.
     * By default, quotes are not handled by a CSVReader.
     * If an empty string is passed to this method, quotes are not handled.
     * @param comment_start indicates that the line should not be processed.
     */
    void
    set_quote(
        char quote
    );

};


}
}

#endif
