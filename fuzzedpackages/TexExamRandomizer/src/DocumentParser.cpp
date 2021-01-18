# include <Rcpp.h>
# include <iostream>
# include <string>
# include <regex>
using namespace Rcpp;
/*
# Author Alejandro Gonzalez Recuenco
# e-mail <alejandrogonzalezrecuenco@gmail.com>
# (C) 2017
*/


class LayerCompiler{
public:
    // const StringVector& x;
    std::vector<std::string> sectionName, commandName;
    std::vector<std::regex> start_regex, end_regex, cmd_regex;
    int last_layer;

    LayerCompiler(
        // const StringVector& _x,
        const StringVector& _sectionName,
        const StringVector& _commandName)
        //: x(_x)
        {

        last_layer = _sectionName.size() - 1;
        for(int i = 0; i< _sectionName.size(); i++){
            sectionName.push_back(as<std::string>(_sectionName[i]));
            commandName.push_back(as<std::string>(_commandName[i]));
            start_regex.push_back(std::regex("^[^%]*\\\\begin\\{"   + sectionName[i] + "\\}"));
            end_regex.push_back  (std::regex("^[^%]*\\\\end\\{"     + sectionName[i] + "\\}"));
            cmd_regex.push_back  (std::regex("^[^%]*\\\\"           + commandName[i] + "([^a-zA-Z]|$)"));

        }

    }

    List compile(const StringVector::const_iterator it_beg, const StringVector::const_iterator it_end, int layer){
        IntegerVector attr_section, attr_command;
        List list;
        int section (0), command (0) ; //command name
        bool isInside (false);
        StringVector::const_iterator it_last (it_beg);



        for(StringVector::const_iterator it_check = it_beg; it_check != it_end; it_check++){

            std::string line_text = as<std::string>(*it_check);

            if (!isInside && isBeginLine(line_text, layer)) {
                if (section == 0){

                    list[namePrior(layer, section, command, attr_section, attr_command)] =
                        range_it(it_last, it_check);
                    // If the first line is the it_last, then this is going to pick both.
                    if (it_last != it_check) {
                        it_last = it_check;
                    } else {
                        it_last = it_check + 1;
                    }
                }
                section++;
                isInside = true;
            } else if (isInside && isCommandLine(line_text, layer)){ // end check

                if (command == 0) {
                    list[nameBegin(layer, section, command, attr_section, attr_command)] =
                        range_it(it_last, it_check);
                } else {
                    if (layer == last_layer){
                        list[nameCommand(layer, section, command, attr_section, attr_command)] =
                            range_it(it_last, it_check);
                    } else {
                        list[nameCommand(layer, section, command, attr_section, attr_command)] =
                            compile(it_last, it_check, layer + 1);
                    }

                }
                command++;

                if (it_last != it_check) {
                    it_last = it_check;
                } else {
                    it_last = it_check + 1;
                }
            } else if (isInside && isEndLine(line_text, layer)) { //command check

                if (command == 0) {
                    list[nameBegin(layer, section, command, attr_section, attr_command)] =
                        range_it(it_last, it_check);
                } else if (layer == last_layer){
                    list[nameCommand(layer, section, command, attr_section, attr_command)] =
                        range_it(it_last, it_check);
                } else {
                    list[nameCommand(layer, section, command, attr_section, attr_command)] =
                        compile(it_last, it_check, layer + 1);
                }
                // it_last shouldnt' be the same of it_check, otherwise the document was ill-formed
                if(it_last ==it_check){
                    stop("Document contains the \\end in the same line as a previous command, \
                         I can't split those lines, stopping compilation ");
                }

                //after the last command we need to save the end on it's own
                list[nameEnd(layer, section, command, attr_section, attr_command)] =
                    range_it(it_check, it_check);
                it_last = it_check + 1;
                //reset the command number after we set everything
                command = 0;
                isInside = false;
            }

        }

        if  (it_last == it_beg){
            list[namePrior(layer, section, command, attr_section, attr_command)] = range_it(it_last, it_end);
            list[namePost(layer, section, command, attr_section, attr_command)] = StringVector::create();
        } else if (it_last == it_end){
            list[namePost(layer, section, command, attr_section, attr_command)] = StringVector::create();
        } else {
            list[namePost(layer, section, command, attr_section, attr_command)] = range_it(it_last, it_end);
        }

        list.attr("section") = attr_section;
        list.attr("command") = attr_command;
        list.attr("section_original") = attr_section;
        list.attr("command_original") = attr_command;


        return list;
    }

private:

    bool isBeginLine(std::string line_of_text, int layer){
        return std::regex_search(line_of_text, start_regex[layer]);
    }
    bool isCommandLine(std::string line_of_text, int layer){
        return std::regex_search(line_of_text, cmd_regex[layer]);
    }
    bool isEndLine(std::string line_of_text, int layer){
        return std::regex_search(line_of_text, end_regex[layer]);
    }


    StringVector range_it(StringVector::const_iterator it_start, StringVector::const_iterator it_end){
        if(it_start > it_end){
            StringVector empty_vector;
            return empty_vector;
        }
        if (it_start == it_end) {
            StringVector vect(it_start, it_start + 1);
            return vect;
        }

        StringVector vect(it_start, it_end);
        return vect;
    }

    void push_sec_cmd(IntegerVector& sec_v, int section, IntegerVector& cmd_v, int command){
        sec_v.push_back(section);
        cmd_v.push_back(command);
    }

    std::string namePrior(int layer, int section, int command, IntegerVector& sec_v,  IntegerVector& cmd_v){
        push_sec_cmd(sec_v, 0, cmd_v,0);
        return "prior_to_" + sectionName[layer];

    }

    std::string nameBegin(int layer, int section, int command, IntegerVector& sec_v,  IntegerVector& cmd_v){
        push_sec_cmd(sec_v, section, cmd_v,0);
        return std::to_string(section) + "_" + sectionName[layer] + "_begin_" + sectionName[layer];

    }
    std::string nameCommand(int layer, int section, int command, IntegerVector& sec_v,  IntegerVector& cmd_v){
        push_sec_cmd(sec_v, section, cmd_v, command);
        return std::to_string(section) + "_" + sectionName[layer] + "_" + std::to_string(command) + "_" + commandName[layer];
    }

    std::string nameEnd(int layer, int section, int command, IntegerVector& sec_v,  IntegerVector& cmd_v){
        push_sec_cmd(sec_v, section, cmd_v, 0);
        return std::to_string(section) + "_" + sectionName[layer] + "_end_" + sectionName[layer];
    }
    std::string namePost(int layer, int section, int command, IntegerVector& sec_v,  IntegerVector& cmd_v){
        push_sec_cmd(sec_v, 0, cmd_v, 0);
        return "post_to_" + sectionName[layer];
    }




};

//' @title Compile Document
//' @description Function that takes a set of lines, \code{x}, that represent a file or a document. And divides it in subsequent layers, structures as a list, as described on the detail section.
//'
//' It assumes \code{x} is representing a 'LaTeX' file that can be compiled as it is before we make any modifications.
//' @inheritParams DivideFile
//' @param layersNames A character vector, with each element representating the environment name to be searched as \code{cmdName} as describe in \code{\link{FindBegin}} and \code{\link{FindEnd}}
//' @param layersCmd A character vector, with the same length as \code{layersNames}. with each element representing the environment command to be serached as \code{cmdName} as described in \code{\link{FindCommand}}.
//'
//'
//' @details
//' Both \code{layersNames} and \code{layersCmd} must have the same length, since for each index, \code{i}, \code{layersNames[i]} and \code{layersCmd[i]} refer to one layer of the tree structure of the document. Consequent layers must be found inside previous layers.
//'
//' If it finds the structure of the document to not be completed, it will throw an error.
//'
//' @return  It returns a list, with each element having a name. Recreating the tree structure identified by  \code{layersNames} and \code{layersCmd} in the text file \code{x}.
//'
//' It first divides the document into two lists:
//' \describe{
//' \item{preamble}{Contains a character vector identifying everything before the \\begin\{document\}}
//' \item{document}{Contains the tree structure identifying the document}
//' }
//'
//' Now, the naming convention for each layer of the document is as follows. We will use the convention \code{<layerName>},  \code{<layerCmd>}.
//'
//'  Note the convention first, everything that it finds prior to the first environment, it throws it into a character vector that it calls \code{prior_to_<layesName>}.
//'  After the first environment \code{<layerName>} ends, it assumes that everything from that \code{\\end\{<layerName>\}} onwards corresponding to the next environment, and it will throw it to the prior part of that one.
//'  \code{post_to_<layerName>}
//' \describe{
//' \item{\code{prior_to_layersName}}{Includes everything up to the first \code{\\begin\{<layerName>} without including that line}
//' \item{\code{1_<layerName>_begin_<layerName>}}{
//'   Includes the \code{\\begin\{layerName\}} for the 1st section, and everything until it finds the first \code{\\<layerCmd>}}
//' \item{\code{1_<layerName>_1_<layerCmd>}}{
//'   Includes everything from the 1\eqn{^{st}} \code{\\<layerCmd>} until the second \code{\\<layerCmd>}, without including the line in which the second command is found
//' }
//' \item{\code{1_<layerName>_2_<layerCmd>}}{
//'  Same thing... and it keeps going until the last  \code{\\<layerCmd>} is found
//' }
//' \item{\code{1_<layerName>_end_<layerName>}}{
//' It includes the \code{\\end\{<layerName>\}} for the 1st section.
//' }
//' \item{...}{
//' It then repeats the same structure for the next environment, changing the naming convention to start with 2_<...> and so on until it does the last environemt}
//' \item{\code{post_to_<layerName>}}{
//' After the last layer ends with \code{\\end\{layerName\}}, it throws the rest of the lines into this last character vector}
//'
//' }
//'
//' This structure is applied recursively to each \code{i_<layerName>_j_<layerCmd>} of the previous layer to find the structure for the next layer.
//' The result is a tree of lists, with names that identify the whole structure, and the ending node of each branch is always a character vector
//'
//' \strong{IMPORTANT NOTE:} Note that this function only rearranges the lines of the document, it can't split a document between a line. So if you want to make sure something always stays together, put them both in the same line. This is intentional, to force a more clear structure on the document that will be parsed
//'
//' In Summary, the sketch of the tree structure would be:
//' \itemize{
//'   \item preamble
//'   \item Document
//'         \itemize{
//'         \item prior_to_LayerName[1]
//'         \item 1_layerName[1]_begin_layerName[1]
//'         \item 1_layerName[1]_1_layerCmd[1]
//'              \itemize{
//'              \item prior_to_LayerName[2]
//'              \item 1_layerName[2]_begin_layerName[2]
//'              \item 1_layerName[2]_1_layerCmd[2]
//'                   \itemize{
//'                   \item Continues...
//'                   }
//'              \item 1_layerName[2]_2_layerCmd[2]
//'                   \itemize{
//'                   \item Continues...
//'                   }
//'              \item ...
//'              \item  post_to_layerName[2]
//'              }
//'         \item 2_layerName[1]_begin_layerName[1]
//'         \item 2_layerName[1]_1_layerCmd[1]
//'              \itemize{
//'              \item ...
//'              }
//'        \item ...
//'        \item n_layerName[1]_end_layerName[1]
//'        \item  post_to_layerName[1]
//'         }
//'
//' }
//'
//' If a  \code{\\<layerCmd>} is not found inside an environment, everything inside that environment is thrown into the begin_layerName part and instead of the numbered environments, an empty character list is added in the middle, with name \code{empty_<layerCmd>} section.
//'
//' @seealso \link{FindStructure} for more information on the details of how the layers are found.
//' @keywords internal
//' @family Structuring Document
// [[Rcpp::export]]
List CompileDocument(
        const StringVector& x,
        const StringVector& layersNames,
        const StringVector& layersCmd){

    if(layersNames.size() != layersCmd.size() || layersNames.size() == 0){
        stop("layersNames and layersCmd don't have the correct format");
    }

    LayerCompiler compiler(layersNames, layersCmd);

    return compiler.compile(x.begin(), x.end(), 0);

}
