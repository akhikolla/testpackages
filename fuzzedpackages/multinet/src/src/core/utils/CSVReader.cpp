#include "core/utils/CSVReader.hpp"
#include "core/exceptions/OperationNotSupportedException.hpp"
#include "core/exceptions/FileNotFoundException.hpp"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>

namespace uu {
namespace core {

// Code taken from anonymous user on stackoverflow - thanks :)
// This reads a line from a stream, and works for all types of endline (\n, \r, \r\n)
std::istream&
uu_getline(std::istream& is, std::string& t)
{
    t.clear();

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for (;;)
    {
        int c = sb->sbumpc();

        switch (c)
        {
        case '\n':
            return is;

        case '\r':
            if (sb->sgetc() == '\n')
            {
                sb->sbumpc();
            }

            return is;

        case EOF:
            // Also handles the case when the last line has no line ending
            is.setstate(std::ios::eofbit);

            if (t.empty())
            {
                is.setstate(std::ios::failbit);
            }

            return is;

        default:
            t += (char)c;
        }
    }
}

// Code adapted from anonymous user on stackoverflow - thanks :)
// This reads a line from a stream up to a field separator,
// and works for all types of endline (\n, \r, \r\n)
std::istream&
uu_getline(std::istream& is, std::string& t, char sep)
{
    t.clear();

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for (;;)
    {
        int c = sb->sbumpc();

        if (c==sep)
        {
            return is;
        }

        switch (c)
        {
        case '\n':
            return is;

        case '\r':
            if (sb->sgetc() == '\n')
            {
                sb->sbumpc();
            }

            return is;

        case EOF:
            // Also handles the case when the last line has no line ending
            is.setstate(std::ios::eofbit);

            if (t.empty())
            {
                is.setstate(std::ios::failbit);
            }

            return is;

        default:
            t += (char)c;
        }
    }
}

CSVReader::
CSVReader()
{
    current_ = "";
    next_ = "";
    has_next_ = true;
    remove_trailing_spaces_ = false;
    next_row_number_ = 0;
    current_row_number_ = 0;
    field_separator_ = ',';
    comment_start_ = "";
    quote_ = '\0';
}


CSVReader::
~CSVReader() {}

void
CSVReader::
open(
    const std::string& path
)
{
    infile_.open(path.data());

    if (infile_.fail())
    {
        throw FileNotFoundException(path);
    }

    // Our code pre-fetches the next line to be read, so that we also know
    // what has_next() should return.
    // If there are comments or empty lines, these lines are automatically skipped
    // (but they are counted when row_number is updated).
    if (comment_start_ == "")
    {
        do
        {
            (uu_getline(infile_, next_))?has_next_=true:has_next_=false;
            next_row_number_++;
        }
        while (next_=="" && has_next_);
    }

    else
    {
        do
        {
            (uu_getline(infile_, next_))?has_next_=true:has_next_=false;
            next_row_number_++;
        }
        while ((next_=="" || next_.find(comment_start_)==0) && has_next_);
    }
}

void
CSVReader::
close()
{
    infile_.close();
}

bool
CSVReader::
has_next(
) const
{
    return has_next_;
}

void
CSVReader::
trim_fields(
    bool value
)
{
    remove_trailing_spaces_ = value;
}


void
CSVReader::
set_expected_num_fields(
    size_t expected_num_fields
)
{
    expected_num_fields_ = expected_num_fields;
}

size_t
CSVReader::
skipped_lines(
) const
{
    return lines_skipped_;
}

std::vector<std::string>
CSVReader::
get_next()
{
    current_ = next_;
    std::istringstream line(next_);
    std::vector<std::string> record;

    while (line)
    {
        std::string field;

        if (!uu_getline(line, field, field_separator_))
        {
            break;
        }

        // If quotation mode is on (quote_ != '\0'), it may happen that two parts of the row now
        // read as different fields will be merged being inside quotation marks.
        // Therefore, if quotation mode is on then the removal of trailing spaces is done
        // later while determining the final fields.
        if (remove_trailing_spaces_ && quote_ == '\0')
        {
            field.erase(field.find_last_not_of(" \t")+1);
            field.erase(0,field.find_first_not_of(" \t"));
        }

        record.push_back(field);
    }

    if (next_.size()>0 && next_.back()==field_separator_)
    {
        record.push_back("");
    }

    // Prepare next line for next call.
    current_row_number_ = next_row_number_;

    // This is the same process done when the file is first opened: we look for
    // the first "valid" row: comments and empty lines are skipped.
    if (comment_start_ == "")
    {
        do
        {
            (uu_getline(infile_, next_))?has_next_=true:has_next_=false;
            next_row_number_++;
        }
        while (next_=="" && has_next_);
    }

    else
    {
        do
        {
            (uu_getline(infile_, next_))?has_next_=true:has_next_=false;
            next_row_number_++;
        }
        while ((next_=="" || next_.find(comment_start_)==0) && has_next_);
    }

    // If quotation mode is on, we may have to aggregate fields that we read as separate
    // fields only because at the previous step quotation marks were not considered, and we
    // may have to processed escaped quotation marks, e.g. by turning a "" (escaped " mark)
    // in the file into a single " in the output field.
    if (quote_ != '\0')
    {
        std::vector<std::string> quoted_record;
        bool quote_on = false;
        size_t quote_start = 0;

        for (size_t pos=0; pos<record.size(); pos++)
        {
            std::string f = record.at(pos);

            if (!quote_on)
            {
                if (f=="" || f.front()!=quote_)
                {
                    if (remove_trailing_spaces_)
                    {
                        f.erase(f.find_last_not_of(" \t")+1);
                        f.erase(0,f.find_first_not_of(" \t"));
                    }

                    quoted_record.push_back(f);
                }

                else
                {
                    quote_on=true;
                    quote_start = pos;
                }
            }

            if (quote_on)
            {
                std::string::size_type n = f.find_last_not_of(quote_);

                if (f!="" && (
                            (n==std::string::npos && quote_start==pos && (f.size())%2==0)
                            ||
                            (n==std::string::npos && quote_start!=pos && (f.size()+1)%2==0)
                            ||
                            (n != std::string::npos && (f.size()-n)%2==0))
                   )
                {
                    quote_on=false;

                    // prepare quoted field

                    std::stringstream ss;

                    for (size_t i=quote_start; i<=pos; i++)
                    {
                        ss << record.at(i);

                        if (i!=pos)
                        {
                            ss << field_separator_;
                        }
                    }

                    std::string composite_field = ss.str().substr(1,ss.str().size()-2);
                    // fixing escaped quotes
                    std::string::size_type n = 0;

                    while ( (n = composite_field.find(escaped_quote_, n)) != std::string::npos )
                    {
                        composite_field.replace(n, 2, quote_as_string_);
                        n += 1;
                    }

                    if (remove_trailing_spaces_)
                    {
                        composite_field.erase(composite_field.find_last_not_of(" \t")+1);
                        composite_field.erase(0,composite_field.find_first_not_of(" \t"));
                    }

                    quoted_record.push_back(composite_field);
                }
            }
        }

        // if an expected number of fields has been specified, this line is skipped and a
        // new one is returned.
        if (expected_num_fields_ != 0 && expected_num_fields_ != quoted_record.size())
        {
            lines_skipped_++;
            return get_next();
        }

        else
        {
            return quoted_record;
        }
    }

    else
    {
        if (expected_num_fields_ != 0 && expected_num_fields_ != record.size())
        {
            lines_skipped_++;
            return get_next();
        }

        else
        {
            return record;
        }
    }
}


std::string
CSVReader::
get_next_raw_line()
{
    current_ = next_;

    // prepare next line for next call
    current_row_number_ = next_row_number_;

    if (comment_start_ == "")
    {
        do
        {
            (uu_getline(infile_, next_))?has_next_=true:has_next_=false;
            next_row_number_++;
        }
        while (next_=="" && has_next_);
    }

    else
    {
        do
        {
            (uu_getline(infile_, next_))?has_next_=true:has_next_=false;
            next_row_number_++;
        }
        while ((next_=="" || next_.find(comment_start_)==0) && has_next_);
    }

    // current_ contains the line that was the next when this method was called
    return current_;
}

std::string
CSVReader::
get_current_raw_line() const
{
    if (current_row_number_==0)
    {
        throw OperationNotSupportedException("cannot retrieve a line from the file before calling get_next()");
    }

    return current_;
}

int
CSVReader::
row_num() const
{
    return current_row_number_;
}

void
CSVReader::
set_field_separator(char separator)
{
    field_separator_ = separator;
}

void
CSVReader::
set_comment(const std::string& comment_start)
{
    comment_start_ = comment_start;
}

void
CSVReader::
set_quote(char quote)
{
    quote_ = quote;
    std::stringstream ss;
    ss << quote_;
    quote_as_string_ = ss.str();
    ss << quote_;
    escaped_quote_ = ss.str();
}

} // namespace core
} // namespace uu
