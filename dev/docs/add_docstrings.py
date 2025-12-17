#!/usr/bin/env python3
# This script was generated with the help of AI (GitHub Copilot)
"""
Add sphinxcontrib-matlab domain compatible docstrings to .m files.
Docstrings are bare-bones (variable names only) for uncertain parameters.
"""

import os
import glob
from pathlib import Path
from collections import Counter


# Greek letter mapping for mathematical rendering in Sphinx
GREEK_LETTERS = {
    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
    'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'omicron', 'pi',
    'rho', 'sigma', 'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega',
}


# Known parameter descriptions
KNOWN_PARAMS = {
    'p': 'Parameter struct',
    'b': 'Design variable bounds struct',
    'X': 'Design variable vector',
    'x': 'Design variable vector',
    'x0': 'Initial design variable vector',
    'X_opt': 'Optimal design variables',
    'X_start_struct': 'Starting design variable struct',
    'obj': 'Objective function value',
    'J': 'Objective function value',
    'flag': 'Optimization output flag',
    'flags': 'Optimization output flags',
    'w': 'Angular wave frequency (rad/s)',
    'k': 'Wavenumber (1/m)',
    'LCOE': 'Levelized cost of energy ($/kWh)',
    'intermed_result_struct': 'Intermediate results struct (cached heavy analyses)',
    'figs': 'Figure handles',
    'fig': 'Figure handle',
    'P_matrix_elec': 'Electrical power matrix',
    'g': 'Constraint vector',
    'val': 'Simulation result values struct',
    'D_d': 'Diameter of damping plate (m)',
    'T_s': 'Draft of spar (m)',
    'failed': 'Indices or names of violated constraints',
    'rho': 'Density of seawater (kg/m^3)',
}


def extract_design_vars_from_bounds(var_bounds_path):
    """Extract design variable names and descriptions from var_bounds.m.
    
    Returns a dict mapping var_name -> description
    """
    var_descriptions = {}
    
    if not os.path.exists(var_bounds_path):
        return var_descriptions
    
    try:
        with open(var_bounds_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Extract var_names array: b.var_names = {...}
        names_start = content.find('b.var_names')
        if names_start == -1:
            return var_descriptions
        
        eq_pos = content.find('=', names_start)
        open_brace = content.find('{', eq_pos)
        close_brace = content.find('}', open_brace)
        names_str = content[open_brace + 1:close_brace]
        
        # Extract var_descs array: b.var_descs = {...}
        descs_start = content.find('b.var_descs')
        if descs_start == -1:
            return var_descriptions
        
        eq_pos = content.find('=', descs_start)
        open_brace = content.find('{', eq_pos)
        close_brace = content.find('}', open_brace)
        descs_str = content[open_brace + 1:close_brace]
        
        # Parse quoted strings from both arrays
        def extract_quoted_strings(text):
            """Extract all single-quoted strings from text."""
            strings = []
            i = 0
            while i < len(text):
                if text[i] == "'":
                    # Find closing quote
                    j = i + 1
                    while j < len(text) and text[j] != "'":
                        j += 1
                    if j < len(text):
                        strings.append(text[i + 1:j])
                    i = j + 1
                else:
                    i += 1
            return strings
        
        names = extract_quoted_strings(names_str)
        descs = extract_quoted_strings(descs_str)
        
        # Map names to descs
        for name, desc in zip(names, descs):
            var_descriptions[name] = desc
    
    except Exception as e:
        print(f"Warning: Could not extract design variables from var_bounds.m: {e}")
    
    return var_descriptions


def extract_parameters_from_parameters_m(parameters_path):
    """Extract parameter names and descriptions from parameters.m table.
    
    Returns a dict mapping param_name -> description
    """
    param_descriptions = {}
    
    if not os.path.exists(parameters_path):
        return param_descriptions
    
    try:
        with open(parameters_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Find all table() calls
        lines = content.split('\n')
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            if line.startswith('table('):
                # Collect the full table() call (might span multiple lines with ...)
                full_call = line
                paren_depth = line.count('(') - line.count(')')
                i += 1
                while paren_depth > 0 and i < len(lines):
                    next_line = lines[i].strip()
                    # Remove MATLAB continuation marker before joining
                    if next_line.endswith('...'):
                        next_line = next_line[:-3].strip()
                    full_call += ' ' + next_line
                    paren_depth += next_line.count('(') - next_line.count(')')
                    i += 1
                
                # Parse table() call to extract name (1st quoted string) and description (6th quoted string)
                # Table structure: table("name", "pretty", {...}, "subsystem", bool, "description", {...})
                # We need to extract quoted strings while respecting bracket depth
                
                quoted_strings = []
                j = 0
                brace_depth = 0
                paren_depth_local = 0
                
                while j < len(full_call):
                    char = full_call[j]
                    
                    # Track bracket depth to skip quoted strings inside {...}
                    if char == '{':
                        brace_depth += 1
                    elif char == '}':
                        brace_depth -= 1
                    elif char == '(':
                        paren_depth_local += 1
                    elif char == ')':
                        paren_depth_local -= 1
                    
                    # Only extract quoted strings that are NOT inside braces or function calls
                    if char == '"' and brace_depth == 0 and paren_depth_local == 1:  # paren_depth_local==1 means inside table()
                        k = j + 1
                        while k < len(full_call) and full_call[k] != '"':
                            k += 1
                        if k < len(full_call):
                            quoted_strings.append(full_call[j + 1:k])
                        j = k + 1
                    else:
                        j += 1
                
                # Table structure yields quoted strings in order:
                # [0]=name, [1]=pretty_name, [2]=subsystem, [3]=description
                if len(quoted_strings) >= 4:
                    param_name = quoted_strings[0]
                    description = quoted_strings[3]
                    param_descriptions[param_name] = description
            else:
                i += 1
    
    except Exception as e:
        print(f"Warning: Could not extract parameters from parameters.m: {e}")
    
    return param_descriptions

# MATLAB docstring template with sphinxcontrib-matlab directives
DOCSTRING_TEMPLATE = '''% {description}
%
% :param {param_name}: {param_desc}
% :returns: {return_desc}'''

DOCSTRING_SINGLE_PARAM = '''% {description}
%
% :param {param_name}: {param_desc}'''

DOCSTRING_NO_PARAM = '''% {description}'''


def extract_function_signature(content):
    """Extract function definition, inputs, and outputs (handles multi-line signatures).
    
    Returns:
        tuple: (func_line, start_idx, end_idx) where end_idx is the line after the signature ends
    """
    lines = content.split('\n')
    func_start_idx = None
    
    # Find the line with 'function' keyword
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if stripped.startswith('function'):
            func_start_idx = i
            break
    
    if func_start_idx is None:
        return None, None, None
    
    # Strategy: Count parentheses () to find end of function inputs
    # Ignore square brackets [] - they're for outputs
    
    func_lines = []
    func_end_idx = func_start_idx
    found_func_open_paren = False  # Track if we've found the function's opening (
    paren_depth = 0
    
    for i in range(func_start_idx, min(func_start_idx + 50, len(lines))):  # max 50 lines for signature
        line = lines[i]
        # Remove inline comments
        line_clean = line.split('%')[0]
        
        # Check for MATLAB continuation marker
        has_continuation = line_clean.rstrip().endswith('...')
        if has_continuation:
            line_clean = line_clean.rstrip()[:-3].rstrip()  # Remove the ...
        
        func_lines.append(line_clean)
        
        # Count parentheses (function inputs) - skip square brackets (outputs)
        for char in line_clean:
            if char == '(':
                paren_depth += 1
                found_func_open_paren = True
            elif char == ')':
                paren_depth -= 1
        
        # If we've found opening paren, counted all parens, and no continuation, we're done
        if found_func_open_paren and paren_depth == 0 and not has_continuation:
            func_end_idx = i
            break
    
    # Join all lines and clean up
    func_line = ' '.join(func_lines).strip()
    # Normalize whitespace
    func_line = ' '.join(func_line.split())
    
    return func_line, func_start_idx, func_end_idx


def is_classdef(content):
    """Check if content is a MATLAB class definition (classdef)."""
    lines = content.split('\n')
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith('classdef '):
            return True
        if stripped.startswith('function '):
            # If we hit function before classdef, it's not a class
            return False
    return False


def extract_classdef_name(content):
    """Extract class name from classdef statement."""
    lines = content.split('\n')
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith('classdef '):
            # Extract name: classdef ClassName
            parts = stripped.split()
            if len(parts) >= 2:
                return parts[1].split('<')[0].split('&')[0].strip()
    return None


def add_classdef_docstring(filepath, design_vars=None, param_descriptions=None):
    """Add or update docstring to a classdef file (class definition).
    
    Adds a marker and docstring right after 'classdef ClassName' line.
    Skips if class already has a docstring.
    """
    if design_vars is None:
        design_vars = {}
    if param_descriptions is None:
        param_descriptions = {}
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        lines = content.split('\n')
        classdef_idx = None
        
        # Find classdef line
        for i, line in enumerate(lines):
            stripped = line.lstrip()
            if stripped.startswith('classdef '):
                classdef_idx = i
                break
        
        if classdef_idx is None:
            return False, "No classdef found"
        
        # Check if already has a docstring after classdef
        # Look for any comment lines (%) or autogenerated marker in next few lines
        has_existing = False
        for j in range(classdef_idx + 1, min(classdef_idx + 10, len(lines))):
            next_line = lines[j].strip()
            
            # If we hit a non-comment, non-whitespace line (like properties), stop
            if next_line and not next_line.startswith('%'):
                break
            
            # If we find comment lines, assume it's a docstring
            if next_line.startswith('%'):
                has_existing = True
                break
        
        if has_existing:
            return False, "Already has docstring"
        
        # Extract class name
        class_name = extract_classdef_name('\n'.join(lines))
        if not class_name:
            class_name = "Class"
        
        # Create simple docstring for class with single marker before
        marker_line = "%[AUTOGENERATED]"
        docstring = f"% {class_name} class definition"
        
        # Insert marker and docstring after classdef line
        lines.insert(classdef_idx + 1, marker_line)
        lines.insert(classdef_idx + 2, docstring)
        
        new_content = '\n'.join(lines)
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
        
        return True, f"Added class docstring for {class_name}"
    
    except Exception as e:
        return False, str(e)



def parse_function_line(func_line):
    """Parse function line to get name, inputs, outputs.
    
    Handles patterns like:
      - function [out1, out2] = func_name(in1, in2)
      - function out = func_name(in1, in2)
      - function func_name(in1, in2)
    Also handles multi-line signatures with ... continuation.
    """
    
    # Remove 'function' keyword and strip whitespace
    if not func_line.startswith('function'):
        return None, None, None, None
    
    remainder = func_line[8:].strip()  # Skip 'function'
    
    # Check if there's an assignment (= sign) for outputs
    if '=' in remainder:
        # Split on = to separate outputs from name+inputs
        output_part, name_and_inputs = remainder.split('=', 1)
        output_part = output_part.strip()
        name_and_inputs = name_and_inputs.strip()
        
        # Parse outputs - could be [out1, out2] or just out
        if output_part.startswith('[') and output_part.endswith(']'):
            outputs_str = output_part[1:-1]  # Remove brackets
            # Remove MATLAB continuation characters (...) from outputs
            outputs_str = outputs_str.replace('...', ' ')
            outputs = [x.strip() for x in outputs_str.split(',') if x.strip()]
        else:
            outputs = [output_part]
    else:
        # No assignment - no output specification
        name_and_inputs = remainder
        outputs = []
    
    # Parse function name and inputs from "name(inputs)"
    if '(' not in name_and_inputs or ')' not in name_and_inputs:
        return None, None, None, None
    
    # Find function name (before opening parenthesis)
    paren_idx = name_and_inputs.index('(')
    func_name = name_and_inputs[:paren_idx].strip()
    
    # Find inputs (between parentheses)
    close_paren_idx = name_and_inputs.rindex(')')
    inputs_str = name_and_inputs[paren_idx + 1:close_paren_idx]
    
    # Remove MATLAB continuation characters (...) from inputs
    # Replace ... with space, then clean up the parameter list
    inputs_str = inputs_str.replace('...', ' ')
    
    # Split by comma and clean up each parameter
    inputs = []
    for param in inputs_str.split(','):
        param = param.strip()
        if param:  # Only add non-empty parameters
            inputs.append(param)
    
    return func_name, inputs, outputs, func_line


def get_param_desc(name, param_descriptions, design_vars):
    """Resolve parameter description with fallbacks.
    Order: parameter table -> design vars -> KNOWN_PARAMS -> greek letter -> pattern xx_over_yy -> name.
    """
    if name in param_descriptions:
        return param_descriptions[name]
    if name in design_vars:
        return design_vars[name]
    if name in KNOWN_PARAMS:
        return KNOWN_PARAMS[name]
    
    # Check for xx_over_yy pattern first (handle Greek letters in both numerator and denominator)
    if '_over_' in name:
        parts = name.split('_over_')
        if len(parts) == 2 and parts[0] and parts[1]:
            # Process numerator for Greek letters
            num_parts = parts[0].split('_')
            if num_parts[0] in GREEK_LETTERS:
                greek_num = f"$\\\\{num_parts[0]}$"
                rest_num = '_'.join(num_parts[1:]) if len(num_parts) > 1 else ''
                numerator = f"{greek_num} {rest_num}".strip()
            else:
                numerator = parts[0]
            
            # Process denominator for Greek letters
            denom_parts = parts[1].split('_')
            if denom_parts[0] in GREEK_LETTERS:
                greek_denom = f"$\\\\{denom_parts[0]}$"
                rest_denom = '_'.join(denom_parts[1:]) if len(denom_parts) > 1 else ''
                denominator = f"{greek_denom} {rest_denom}".strip()
            else:
                denominator = parts[1]
            
            return f"{numerator} / {denominator}".replace('_', ' ')
    
    # Check for Greek letters at start (e.g., theta, gamma_phase_f)
    parts = name.split('_')
    if parts[0] in GREEK_LETTERS:
        greek_part = parts[0]
        rest = '_'.join(parts[1:]) if len(parts) > 1 else ''
        greek_symbol = f"$\\\\{greek_part}$"
        if rest:
            return f"{greek_symbol} {rest}".replace('_', ' ')
        return greek_symbol
    
    return name



def generate_docstring(func_name, inputs, outputs, design_vars=None, param_descriptions=None):
    """Generate docstring based on function signature.
    
    Args:
        func_name: Name of the function
        inputs: List of input parameter names
        outputs: List of output parameter names
        design_vars: Dict mapping design var names to their descriptions
        param_descriptions: Dict mapping parameter names to descriptions from parameters.m
    
    Returns the docstring text only (marker goes before function, not in docstring).
    """
    if design_vars is None:
        design_vars = {}
    if param_descriptions is None:
        param_descriptions = {}
    
    description = f"Function {func_name}"
    
    if len(outputs) == 0:
        # Just description, no returns
        if len(inputs) == 0:
            return f"% {description}"
        else:
            param_lines = []
            for param in inputs:
                # Check sources in priority order: param_descriptions, design_vars, KNOWN_PARAMS
                desc = get_param_desc(param, param_descriptions, design_vars)
                param_lines.append(f"% :param {param}: {desc}")
            return f"% {description}\n%\n" + '\n'.join(param_lines)
    
    # Build full docstring with inputs and outputs
    docstring = f"% {description}\n%"
    
    # Add input parameters
    for param in inputs:
        desc = get_param_desc(param, param_descriptions, design_vars)
        docstring += f"\n% :param {param}: {desc}"
    
    # Add return values
    for ret in outputs:
        desc = get_param_desc(ret, param_descriptions, design_vars)
        docstring += f"\n% :returns: {desc}"
    
    return docstring



def has_docstring(content):
    """Check if function has docstring, return (has_manual, has_autogenerated).
    
    Returns:
        tuple: (has_manual_docstring, has_autogenerated_docstring)
    
    A docstring is considered "manual" only if it exists WITHOUT an autogenerated marker.
    If there's an autogenerated marker, it should always be replaced.
    """
    lines = content.split('\n')
    for i, line in enumerate(lines):
        stripped = line.lstrip()
        if stripped.startswith('function'):
            # Check for marker on line before function
            has_autogen_marker = False
            if i > 0 and '%[AUTOGENERATED]' in lines[i-1]:
                has_autogen_marker = True
            
            # If autogenerated marker exists, always return (False, True) to indicate
            # we should replace the docstring
            if has_autogen_marker:
                return (False, True)
            
            # Otherwise, check if there's any docstring at all
            has_manual = False
            # Only treat as manual if sphinxcontrib-matlab directives present after function
            for j in range(i + 1, min(i + 40, len(lines))):
                next_line = lines[j].strip()
                if ':param' in next_line or ':returns' in next_line:
                    has_manual = True
                    return (has_manual, False)
                if next_line and not next_line.startswith('%'):
                    return (has_manual, False)

            return (has_manual, False)
    return (False, False)



def remove_autogenerated_docstring(content):
    """Remove only the autogenerated docstring lines after function.
    
    Removes ALL occurrences of autogenerated docstring elements, including duplicates:
    - Lines with % :param
    - Lines with % :returns
    - Lines with % Function
    - ONLY the bare % lines that are part of autogenerated docstring structure
    
    Preserves:
    - Any other comment content (manual documentation)
    - Bare % lines within manual documentation
    - Code
    - Empty lines that are not part of docstring
    """
    lines = content.split('\n')
    result = []
    found_function = False
    in_autogen_docstring = False
    last_line_was_autogen = False
    i = 0
    
    while i < len(lines):
        line = lines[i]
        stripped = line.lstrip()
        
        # Skip all %[AUTOGENERATED] markers
        if '%[AUTOGENERATED]' in stripped:
            i += 1
            continue
        
        # If we found a function line, mark it and handle signature (may span multiple lines with ...)
        if stripped.startswith('function '):
            result.append(line)
            found_function = True
            in_autogen_docstring = False
            last_line_was_autogen = False
            i += 1
            
            # Handle multi-line function signatures with ... continuation
            while i < len(lines):
                next_line = lines[i]
                next_stripped = next_line.lstrip()
                
                # If it's a comment (starts with %), we've reached the docstring, exit loop
                if next_stripped.startswith('%'):
                    break
                
                # If it's empty, include it and continue
                if not next_stripped:
                    result.append(next_line)
                    i += 1
                    continue
                
                # Remove inline comments to check for continuation marker
                line_without_comment = next_line.split('%')[0].rstrip()
                
                # If line ends with continuation marker, it's part of the signature
                if line_without_comment.endswith('...'):
                    result.append(next_line)
                    i += 1
                else:
                    # No continuation marker - this is the last line of the signature
                    result.append(next_line)
                    i += 1
                    # Now we're done with signature, exit loop to process docstring
                    break
            
            continue
        
        # After function signature ends, selectively remove autogenerated docstring lines
        if found_function:
            # Stop processing if we hit non-comment, non-empty content
            if stripped and not stripped.startswith('%'):
                found_function = False
                in_autogen_docstring = False
                last_line_was_autogen = False
                result.append(line)
                i += 1
                continue
            
            # Empty lines - only skip if they're between autogenerated docstring elements
            if not stripped:
                # Only skip if we are currently in the docstring AND next line is also docstring
                if in_autogen_docstring and i + 1 < len(lines):
                    next_stripped = lines[i + 1].lstrip()
                    # Check if next line is part of autogenerated docstring
                    next_is_autogen = (
                        next_stripped.startswith('% :param') or 
                        next_stripped.startswith('% :returns') or
                        next_stripped.startswith('% Function ') or
                        next_stripped == '%'
                    )
                    if next_is_autogen:
                        # Skip this empty line (it's between docstring elements)
                        i += 1
                        continue
                
                # Keep empty lines that are not between docstring elements
                result.append(line)
                last_line_was_autogen = False
                i += 1
                continue
            
            # Identify autogenerated docstring elements
            is_function_line = stripped.startswith('% Function ')
            is_param_line = stripped.startswith('% :param')
            is_returns_line = stripped.startswith('% :returns')
            is_bare_percent = stripped == '%'
            
            # A bare % is only autogenerated if it's between Function/params/returns
            skip_bare_percent = False
            if is_bare_percent and in_autogen_docstring:
                # Look ahead to see if next is still autogenerated content
                if i + 1 < len(lines):
                    next_stripped = lines[i + 1].lstrip()
                    next_is_autogen = (
                        next_stripped.startswith('% :param') or 
                        next_stripped.startswith('% :returns') or
                        next_stripped.startswith('% Function ') or
                        next_stripped == '%'
                    )
                    if next_is_autogen:
                        skip_bare_percent = True
                    else:
                        # End of autogenerated docstring
                        in_autogen_docstring = False
            
            # Track when we're in autogenerated docstring section
            if is_function_line or is_param_line or is_returns_line:
                in_autogen_docstring = True
            
            # Decide whether to skip this line
            is_autogen_line = (
                is_function_line or 
                is_param_line or 
                is_returns_line or 
                skip_bare_percent
            )
            
            if is_autogen_line:
                # Skip autogenerated lines
                last_line_was_autogen = True
                i += 1
                continue
            else:
                # Keep everything else (manual documentation)
                result.append(line)
                last_line_was_autogen = False
                i += 1
                continue
        
        # Before function line, keep everything
        result.append(line)
        last_line_was_autogen = False
        i += 1
    
    return '\n'.join(result)


def count_param_usage(filepaths, exclusions=None):
    """Count parameter name occurrences (inputs + outputs) across MATLAB files.

    exclusions: set of names to skip (e.g., known params, design vars, parameter table, '~').
    """
    if exclusions is None:
        exclusions = set()
    counter = Counter()
    for filepath in filepaths:
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            func_line, _, _ = extract_function_signature(content)
            if func_line is None:
                continue
            func_name, inputs, outputs, _ = parse_function_line(func_line)
            if func_name is None:
                continue
            for name in inputs + outputs:
                if name and name not in exclusions:
                    counter[name] += 1
        except Exception:
            # Ignore files we cannot parse for counting purposes
            continue
    return counter


def add_docstring_to_file(filepath, design_vars=None, param_descriptions=None):
    """Add or update docstring to a .m file.
    
    Replaces autogenerated docstrings, skips files with manual docstrings.
    Places single marker before function line.
    """
    if design_vars is None:
        design_vars = {}
    if param_descriptions is None:
        param_descriptions = {}

    # Explicit skip list for files with long manual headers we should not touch
    skip_files = {
        str(Path('mdocean') / 'plots' / 'util' / 'table2latex.m')
    }
    if any(filepath.endswith(s) for s in skip_files):
        return False, "Already has manual docstring"
    
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Check if already has docstring
        has_manual, has_autogen = has_docstring(content)
        
        if has_manual:
            return False, "Already has manual docstring"
        
        # Extract function signature BEFORE modifying content
        func_line, line_start_idx, line_end_idx = extract_function_signature(content)
        if func_line is None:
            return False, "Could not parse function"
        
        # Now remove the old autogenerated docstring if present
        if has_autogen:
            # Remove the old autogenerated docstring
            content = remove_autogenerated_docstring(content)
            # Re-extract signature after removal (indices may have shifted)
            func_line, line_start_idx, line_end_idx = extract_function_signature(content)
        
        # Parse function signature
        func_name, inputs, outputs, _ = parse_function_line(func_line)
        if func_name is None:
            return False, "Could not parse function signature"
        
        # Generate docstring (just the comment block, marker goes separately)
        docstring = generate_docstring(func_name, inputs, outputs, design_vars, param_descriptions)
        
        # Insert marker before function, docstring after function signature
        lines = content.split('\n')
        lines.insert(line_start_idx, '%[AUTOGENERATED]')
        # After insertion, indices shift by 1
        lines.insert(line_end_idx + 2, docstring)
        
        new_content = '\n'.join(lines)
        
        # Write back
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_content)
        
        return True, f"Added docstring (inputs: {inputs}, outputs: {outputs})"
    
    except Exception as e:
        return False, str(e)


def main():
    """Process all .m files in mdocean folder."""
    # Resolve mdocean path relative to this script (dev/docs/add_docstrings.py)
    script_dir = Path(__file__).parent.parent.parent  # Go up from dev/docs/
    mdocean_path = script_dir / 'mdocean'
    
    # Extract design variable descriptions from var_bounds.m
    var_bounds_file = mdocean_path / 'inputs' / 'var_bounds.m'
    design_vars = extract_design_vars_from_bounds(str(var_bounds_file))
    print(f"Loaded {len(design_vars)} design variable descriptions")
    
    # Extract parameter descriptions from parameters.m
    parameters_file = mdocean_path / 'inputs' / 'parameters.m'
    param_descriptions = extract_parameters_from_parameters_m(str(parameters_file))
    print(f"Loaded {len(param_descriptions)} parameter descriptions\n")
    
    # Find all .m files, excluding generated, SAFE, and OpenFLASH folders
    pattern = str(mdocean_path / '**' / '*.m')
    m_files = glob.glob(pattern, recursive=True)
    
    # Filter out generated, SAFE, and OpenFLASH
    m_files = [f for f in m_files if '/generated/' not in f and '/SAFE/' not in f and '/OpenFLASH/' not in f]
    
    print(f"Found {len(m_files)} .m files to process\n")
    
    success_count = 0
    skip_count = 0
    error_count = 0
    
    for i, filepath in enumerate(sorted(m_files), 1):
        rel_path = os.path.relpath(filepath, mdocean_path)
        
        # Check if it's a classdef file
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            is_class = is_classdef(content)
        except:
            is_class = False
        
        # Handle classdef or regular function
        if is_class:
            success, msg = add_classdef_docstring(filepath, design_vars, param_descriptions)
        else:
            success, msg = add_docstring_to_file(filepath, design_vars, param_descriptions)
        
        if success:
            success_count += 1
            print(f"[{i}/{len(m_files)}] ✓ {rel_path}")
            print(f"         {msg}\n")
        else:
            if "manual docstring" in msg:
                skip_count += 1
            else:
                error_count += 1
            if i % 20 == 0 or error_count > 0:  # Print progress periodically
                print(f"[{i}/{len(m_files)}] - {rel_path}: {msg}")
    
    print("\n" + "="*60)
    print(f"Summary:")
    print(f"  Added docstrings: {success_count}")
    print(f"  Already had docstrings: {skip_count}")
    print(f"  Errors: {error_count}")
    print(f"  Total processed: {len(m_files)}")

    # Show most common parameter names to inform KNOWN_PARAMS additions
    exclusion_set = set(KNOWN_PARAMS.keys()) | set(design_vars.keys()) | set(param_descriptions.keys()) | {'~'}
    usage = count_param_usage(m_files, exclusions=exclusion_set)
    if usage:
        print("\nTop parameter names (inputs/outputs):")
        for name, count in usage.most_common(20):
            print(f"  {name}: {count}")



if __name__ == '__main__':
    main()
