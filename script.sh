version=v09.1
Nelectrons=6
Tds=(0.3 1.0 2.0 3.0 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5 1.7 4.6 8.3)
#Tds=(7.0)
for Td in ${Tds[@]}; do
	echo Release/ThreeDimSimulation 42 ${Td} ${Nelectrons} "Output/${version}/eData_${Td}Td.root"
	Release/ThreeDimSimulation 42 ${Td} ${Nelectrons} "Output/${version}/eData_${Td}Td.root"
done

