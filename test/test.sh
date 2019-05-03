
#!/bin/bash


../convino ../examples/differential_example.txt --prefix diffex1

../convino ../examples/exampleconfig.txt --prefix ex1

../differential_example






echo "------------------------------ diff text based example ---------------------------------------"
diff test_ex1_result.txt ex1_result.txt
echo "-------------------------- diff c++  based differential example ---------------------"
diff differentialExample_output.txt test_differentialExample_output.txt
echo "-------------------------- diff text based differential example (small deviations expected b/c normalisation)------------------------------ "
diff test_diffex1_result.txt diffex1_result.txt