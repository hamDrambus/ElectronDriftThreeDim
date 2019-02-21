version=v06.1
Nelectrons=1
Tds=(7.0)
for Td in ${Tds[@]}; do
	echo Release/ThreeDimSimulation 42 ${Td} ${Nelectrons} "Output/${version}/eData_${Td}Td.root"
	Release/ThreeDimSimulation 42 ${Td} ${Nelectrons} "Output/${version}/eData_${Td}Td.root"
done

