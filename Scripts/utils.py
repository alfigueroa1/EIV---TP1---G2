import numpy as np

def to_unit(num, sf=3):
	"""Round a number to the specified significant figures (default is 3) and add the corresponding unit prefix.

	Args:
		var (float): number to be formatted.
		sf (int, optional): Significant figures. Defaults to 3.

	Returns:
		string: formatted number and unit prefix.
	"""
	units = ['', 'k', 'M', 'G', 'T', 'p', 'n', 'u', 'm']
	e_count = 0
	value = num
	while e_count <= 4 and e_count >= -4:
		if 1 <= value < 1E3:
			break
		if value >= 1E3:
			value /= 1E3
			e_count += 1
		if value < 1:
			value *= 1E3
			e_count -= 1
	
	if e_count > 4 or e_count < -4:
		return '{:.{}g}'.format(num, sf-1).rstrip('0').rstrip('.')
	else:
		return np.format_float_positional(value, precision=sf, unique=True, fractional=False, trim='-') + units[e_count]
	