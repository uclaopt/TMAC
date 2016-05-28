#! /bin/bash

exec_path="./bin/"
data=""
epochs=""
num_threads=""
label=""


# function to display welcome message
function welcome_msg {
    echo "------------------------------------------"
    echo "|        Welcome to ARock demos          |"
    echo "------------------------------------------"
    echo "| Enter 0 to exit                        |"    
    echo "| Enter 1 to run least squares           |"
    echo "| Enter 2 to run lasso                   |"
    echo "| Enter 3 to run l2 svm                  |"
    echo "| Enter 4 to run l1 logistic regression  |"
    echo "| Enter 5 to run demo PRS                |"   
    echo "------------------------------------------"    
}




# function to proess the input from the terminal
function process_app_with_data {
    echo "Please enter the data file name:"
    read data </dev/tty
    echo "Please enter the label file name:"
    read label </dev/tty
    echo "Please enter the number epochs:"
    read epochs </dev/tty
    echo "Please enter the number of threads:"
    read num_threads </dev/tty
}

function process_app_with_no_data {
    echo "please enter number of epochs"
}



#################
welcome_msg

while read input </dev/tty; do
    case $input in
	0 )
	    exit 0 ;;
	1 )
	    app="arock_gd_ls"
	    process_app_with_data
	    exec="$exec_path$app -data ./data/$data -label ./data/$label -epoch $epochs -nthread $num_threads"
	    echo "start ARock to solve the least squares problems."
	    $exec
	    exit 0 ;;
	2 )
	    app="arock_fbs_lasso"
	    process_app_with_data
	    exec="$exec_path$app -data ./data/$data -label ./data/$label -epoch $epochs -nthread $num_threads"
	    echo "starting ARock to solve LASSO"
	    $exec
	    exit 0 ;;
	3 )
	    app="arock_fbs_l2_svm"
	    process_app_with_data
	    exec="$exec_path$app -data ./data/$data -label ./data/$label -epoch $epochs -nthread $num_threads"
	    echo "starting ARock to solve L2_SVM"
	    $exec
	    exit 0 ;;
	4 )
	    app="arock_fbs_l1_log"
	    process_app_with_data
	    exec="$exec_path$app -data ./data/$data -epoch $epochs -nthread $num_threads"
	    echo "starting ARock to solve sparse logistic regression"
	    $exec
	    exit 0 ;;
	5 )
	    app="arock_prs_demo"
	    exec="$exec_path$app"
	    echo "starting ARock to demo PRS"
	    $exec
	    exit 0 ;;
	
    esac

done

