#include<string>
#include<vector>
#include<fstream>
#include<iostream>
#include<sstream>
#include<map>

#include"Vector.h"
#include"Matrix.h"

// identifiers
const std::string operator_process_identifier = "#process";
const std::string create_identifier = "#create";

/*
	[0] = splitting_scheme type
	[1] = splitting_scheme name
	[2] = forward_operator type
	[3] = forward_operator template param
	[4] = forward_operator name
	[5] = backward_operator type
	[6] = backward_operator template param
	[7] = bakward_operator name
*/	
std::vector<std::string> the_holder_of_all = {
	"SPLITTING SCHEME TYPE UNCHANGED",
	"SPLITTING SCHEME NAME UNCHANGED",
	"FORWARD OPERATOR TYPE UNCHANGED",
	"FORWARD OPERATOR TEMPLATE PARAM UNCHANGED",
	"FORWARD OPERATOR NAME UNCHANGED",
	"BACKWARD OPERATOR TYPE UNCHANGED",
	"BACKWARD OPERATOR TEMPLATE PARAM UNCHANGED",
	"BACKWARD OPERATOR NAME UNCHANGED",  
};

// a sort of dictionary of templatized operators
const std::vector<std::string> template_operators = {
	"forward_grad_for_square_loss",
	"forward_grad_for_log_loss",
	"forward_grad_for_square_hinge_loss",
	"forward_grad_for_qp",
	"forward_grad_for_dual_svm",
	"linear_eqn_jacobi_operator",
	"forward_grad_for_huber_loss"
};

// operator map that links operators to their respective templates
std::map<std::string, std::string> operator_map;

// function declarations
std::string process_operator_line(std::string line, bool working_with_forward_operator, std::string operator_type, std::string operator_name);
std::string create_line(std::string line);

int main() 
{
	///////////////////////////// Set up Operator Map /////////////////////////////
	std::map<std::string, std::string> om;
	
	// Group 1 Operators
	om["prox_l1"] = "template_1.txt";
	om["prox_sum_square"] = "template_1.txt";
	om["prox_l2"] = "template_1.txt";
	om["prox_log_barrier"] = "template_1.txt";
	om["proj_positive_cone"] = "template_1.txt";
	om["proj_prob_simplex"] = "template_1.txt";

	// Group 2 Operators
	om["prox_huber"] = "template_2.txt";

	// Group 3 Operators
	om["prox_elastic_net"] = "template_3.txt";

	// Group 4 Operators
	om["proj_box"] = "template_4.txt";

	// Group 5 Operators
	om["proj_l1_ball"] = "template_5.txt";
	om["proj_l2_ball"] = "template_5.txt";

	// Group 6 Operators
	om["proj_hyperplane"] = "template_6.txt";

	// Group 7 Operators
	om["forward_grad_for_square_loss"] = "template_7.txt";
	om["forward_grad_for_log_loss"] = "template_7.txt";
	om["forward_grad_for_square_hinge_loss"] = "template_7.txt";

	// Group 8 Operators
	om["forward_grad_for_qp"] = "template_8.txt";

	// Group 9 Operators
	om["forward_grad_for_dual_svm"] = "template_9.txt";

	// Group 10 Operators
	om["linear_eqn_jacobi_operator"] = "template_10.txt";

	// Group 11 Operators
	om["forward_grad_for_huber_loss"] = "template_11.txt";

	// now set equal to global variable map
	operator_map = om;

	///////////////////////////// Input File Reading /////////////////////////////

	//open the initial template (found in a txt file, identifiers strategically placed)
	std::ifstream in("main_template.txt");
	//create storage for each of the lines in the template (the txt file)
	std::vector<std::string> lines;
	//temp variable for reading from stream
	std::string temp;

	//grab lines until grab lines fails
	while (getline(in, temp)) {
		lines.push_back(temp);
	}

	// we process all the lines for #create
	for (auto& line : lines) {
		//Key assumption, only one #create per line
		size_t loc = line.find(create_identifier);
		if (loc != std::string::npos) { // if we CAN find the identifier...
			line = create_line(line.substr(loc + create_identifier.size()));
		}
	}

	///////////////////////////// Output File Construction /////////////////////////////

	std::cout << "Constructing the following file..." << std::endl;
	
	for (auto& line : lines) 
	{
		std::cout << line << std::endl; // plop everything out to the console
	}
	
	std::ofstream out("C:\\Users\\Elliot\\Documents\\Visual Studio 2015\\Projects\\target\\target\\test_file.cpp");
	// can be replaced with full path of computer to the main of another project
	// specifying the output

	for (auto& line : lines) 
	{
		out << line << std::endl; // now we throw everything into the specified file
	}
	out.close();

	return 1;
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////


std::string process_operator_line(std::string line, bool working_with_forward_operator, std::string operator_type, std::string operator_name) {
	//extracting strings
	std::stringstream iss(line);
	//formatting output
	std::stringstream oss;
	//extract next two string (NOTE: ONLY THOSE NEXT TWO STRINGS)
	std::string type, name;
	//#process| type name
	iss >> type >> name;

	if (type == "operator")
	{
		// need operator_type param here for this initial statement
		std::cout << "Detected " << operator_type << " " << operator_name << " in template. " << std::endl;

		// now switch to local variables instead of func. params.
		type = operator_type;
		name = operator_name;

		// cross check to see if the operator given is a template operator or not.
		for (const auto& op : template_operators)
		{
			if (type == op)
			{
				std::string template_parameter;
				std::cout << "You have selected a template operator. What would you like the template parameter to be? Please input: ";
				std::cin >> template_parameter;

				// this will come in handy when outputting Matrix pointers...
				if (working_with_forward_operator)
					the_holder_of_all[3] = template_parameter;
				else
					the_holder_of_all[6] = template_parameter;
				
				// now we set type to be it's full form w/ template type
				type = type + "<" + template_parameter + ">";

				// now store in either forward_operator_type spot or backward_operator_type spot
				if (working_with_forward_operator)
				{
					the_holder_of_all[2] = type;
				}
				else
				{
					the_holder_of_all[5] = type;
				}
			}
		}

		oss << type << " " << name << ";";
	}
	else if (type == "double") // if double
	{
		std::cout << "Detected " << type << " " << name << " in template. ";

		std::cout << "What value do you want double " << name << " to be? Please input: ";
		double val;
		std::cin >> val;
		oss << operator_name << "." << name << " = " << val << ";";
	}
	else if (type == "Vector") // if Vector
	{
		std::cout << "Detected " << type << " " << name << " in template. ";

		std::cout << "What file do you want to use to construct " << name << "? Please input: ";
		std::string filename;
		std::cin >> filename;
		oss << operator_name << "." << name << ".set_filename(" << filename << ");";
	}
	else if (type == "Matrix") // if Matrix
	{
		std::cout << "Detected " << type << " " << name << " in template. ";

		std::cout << "What file do you want to use to construct " << name << "? Please input: ";
		std::string filename;
		std::cin >> filename;
		oss << operator_name << "." << name << ".set_filename(" << filename << ");";
	}
	else if (type[type.length() - 1] == '*') // if last character is a star
	{
		// now here the thing is that if we have a Mat, we need to replace it...
		if (type == "Mat*")
		{
			if (working_with_forward_operator)
				type = the_holder_of_all[3] + "*";
			else
				type = the_holder_of_all[6] + "*";
		}

		std::cout << "Detected " << type << " " << name << " in template. ";

		std::cout << "What file do you want to use to construct " << name << "? Please input: ";
		std::string filename;
		std::cin >> filename;

		oss << operator_name << "." << name << " = new " << type.substr(0, type.length() - 1) << "(\"" << filename << "\");";
	}
	else
	{
		std::cout << "unrecognized type. probably should exit." << std::endl;
	}
	return oss.str();
}

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

std::string create_line(std::string line) 
{
	std::stringstream iss(line);
	std::stringstream oss;
	std::string tag;

	// #create tag
	iss >> tag;

	if (tag == "TMAC")
	{
		// The usual.
		std::ifstream in("TMAC_template.txt");
		std::vector<std::string> lines;
		std::string temp;

		while (getline(in, temp)) 
		{
			lines.push_back(temp);
		}

		// we process all the lines for #create
		for (auto& line : lines) 
		{
			//Key assumption, only one #create per line
			size_t loc = line.find(create_identifier);

			if (loc != std::string::npos) 
			{
				line = create_line(line.substr(loc + create_identifier.size()));
			}

			oss << line << std::endl;
		}

		oss << "TMAC(" << the_holder_of_all[1] << ");";
	}
	else if (tag == "splitting_scheme")
	{
		// Prompt user for splitting scheme
		std::string scheme, scheme_template;
		std::cout << "Give me a splitting_scheme: ";
		std::cin >> scheme;

		// Scheme selection branching
		if (scheme == "FBS" || scheme == "fbs")
		{
			scheme_template = "fbs_template.txt";
			the_holder_of_all[0] = "ForwardBackwardSplitting";
		}

		// The usual.
		std::ifstream in(scheme_template);
		std::vector<std::string> lines;
		std::string temp;

		while (getline(in, temp)) 
		{
			lines.push_back(temp);
		}

		for (auto& line : lines) 
		{
			size_t loc = line.find(create_identifier);
			
			if (loc != std::string::npos) 
			{
				line = create_line(line.substr(loc + create_identifier.size()));
			}

			oss << line << std::endl;
		}

		// Give generic name to splitting scheme
		the_holder_of_all[1] = "banana_split";

		// QUESTION: Are we assuming the user isn't stupid and didn't name his operators the same as his scheme?
		if (the_holder_of_all[0] == "ForwardBackwardSplitting")
		{
			// change the splitting_scheme name
			oss << the_holder_of_all[0] + "<" + the_holder_of_all[2] + "," + the_holder_of_all[5] + "> " + the_holder_of_all[1] + "("
				+ the_holder_of_all[4] + "," + the_holder_of_all[7] + ");";
		} // will add to this if/else structure as we add additional schemes...

	}
	else if (tag == "forward")
	{
		// Prompt for a foward operator
		std::string forward_operator;
		std::cout << "What kind of forward operator would you like to use: ";
		std::cin >> forward_operator;

		// provide generic name
		std::string forward_operator_name = "fop";

		// Store forward operator type, name and template group
		the_holder_of_all[2] = forward_operator;
		the_holder_of_all[4] = forward_operator_name;
		std::string chosen_template = operator_map[forward_operator];

		// The usual.
		std::ifstream in(chosen_template);
		std::vector<std::string> lines;
		std::string temp;

		while (getline(in, temp)) 
		{
			lines.push_back(temp);
		}

		for (auto& line : lines) 
		{
			size_t loc = line.find(operator_process_identifier);

			if (loc != std::string::npos) 
			{ 
				// NOTE: we PROCESS, not create...
				line = process_operator_line(line.substr(loc + operator_process_identifier.size()), true, forward_operator, forward_operator_name);
			}

			oss << line << std::endl;
		}
	}
	else if (tag == "backward")
	{
		// Prompt user for backward operator
		std::string backward_operator;
		std::cout << "What kind of backward operator would you like to use: ";
		std::cin >> backward_operator;

		// generic name for back operator
		std::string backward_operator_name = "bop";

		// Store backward operator type, name and template group
		the_holder_of_all[5] = backward_operator;
		the_holder_of_all[7] = backward_operator_name;
		std::string chosen_template = operator_map[backward_operator];

		// The usual.
		std::ifstream in(chosen_template);
		std::vector<std::string> lines;
		std::string temp;

		while (getline(in, temp)) 
		{
			lines.push_back(temp);
		}

		for (auto& line : lines)
		{
			size_t loc = line.find(operator_process_identifier);

			if (loc != std::string::npos)
			{
				// NOTE: we PROCESS, not create...
				line = process_operator_line(line.substr(loc + operator_process_identifier.size()), false, backward_operator, backward_operator_name);
			}

			oss << line << std::endl;
		}
	}

	return oss.str();
}
