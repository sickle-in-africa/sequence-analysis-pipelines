#
#	GENERATE JSON - python tool
#
#	* creates a new json file from an old
#	* one, by modifying a value corresponding
#	* to a key supplied by the user. 
#	*
#	* The input and output json file paths can
#	* be the same, in which case the input json
#	* file is over-written.
#	* If a key is given for which there is no
#	* corresponding key in the input json, a new
#	* key value pair will be created.
#
#	Jack Morrice
#
##########################################

import sys, argparse
import json


def main():

	args=get_shell_arguments()

	with open(args.input_json) as json_file:

		# read input json file
		inputs = json.load(json_file)

		# modify the chosen value
		inputs[args.key] = args.new_value

		# write out to file
		output=open(args.output_json, 'w')
		output.write(json.dumps(inputs, indent=2))
		output.close()


#
#  functions
#
def get_shell_arguments():

	parser = argparse.ArgumentParser()

	parser.add_argument(
		'--input_json',
		type=str,
		help='input json file to modify')

	parser.add_argument(
		'--key',
		type=str,
		help='key of value to modify')

	parser.add_argument(
		'--new_value',
		type=str,
		help='new value to replace the input one at <key> location')

	parser.add_argument(
		'--output_json',
		type=str,
		help='output json file name to save to')

	return parser.parse_args()


if __name__ == '__main__':
	main()