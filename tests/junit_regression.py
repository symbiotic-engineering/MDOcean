import xml.etree.ElementTree as ET
import sys
from datetime import datetime
from glob import glob

def parse_junit_xml(file_path):
    tree = ET.parse(glob(file_path))
    root = tree.getroot()
    test_cases = {}
    test_elements = {}  # Store actual XML elements
    for testsuite in root.findall('testsuite'):
        for testcase in testsuite.findall('testcase'):
            name = testcase.get('name')
            classname = testcase.get('classname')
            key = f"{classname}.{name}"
            test_cases[key] = testcase.attrib
            test_elements[key] = testcase
    return test_cases, test_elements, root

def get_test_status(test_element):
    """Determine if a test passed or failed based on child elements."""
    # Check for failure, error, or skipped elements
    if test_element.find('failure') is not None:
        return 'FAILED'
    elif test_element.find('error') is not None:
        return 'ERROR'
    elif test_element.find('skipped') is not None:
        return 'SKIPPED'
    else:
        return 'PASSED'

def create_filtered_junit_xml(base_file, compare_file, output_file):
    base_cases, base_elements, base_root = parse_junit_xml(base_file)
    compare_cases, compare_elements, compare_root = parse_junit_xml(compare_file)
    
    # Create new XML structure
    new_root = ET.Element('testsuites')
    new_testsuite = ET.SubElement(new_root, 'testsuite')
    
    # Copy attributes from original testsuite if available
    original_testsuite = compare_root.find('testsuite')
    if original_testsuite is not None:
        for attr, value in original_testsuite.attrib.items():
            new_testsuite.set(attr, value)
    
    new_or_changed_tests = []
    
    # Find new test cases
    for key in compare_cases:
        if key not in base_cases:
            new_or_changed_tests.append((key, compare_elements[key], "NEW"))
    
    # Find changed test cases (only pass/fail status changes)
    for key in compare_cases:
        if key in base_cases:
            base_status = get_test_status(base_elements[key])
            compare_status = get_test_status(compare_elements[key])
            if base_status != compare_status:
                new_or_changed_tests.append((key, compare_elements[key], "CHANGED"))
    
    # Add filtered test cases to new XML
    for key, test_element, status in new_or_changed_tests:
        # Create a copy of the test element
        new_testcase = ET.SubElement(new_testsuite, 'testcase')
        for attr, value in test_element.attrib.items():
            new_testcase.set(attr, value)
        
        # Copy child elements (like failure, error, system-out, etc.)
        for child in test_element:
            new_testcase.append(child)
        
        # Add a property to indicate the status
        properties = new_testcase.find('properties')
        if properties is None:
            properties = ET.SubElement(new_testcase, 'properties')
        
        property_elem = ET.SubElement(properties, 'property')
        property_elem.set('name', 'regression_status')
        property_elem.set('value', status)
    
    # Update testsuite attributes
    new_testsuite.set('tests', str(len(new_or_changed_tests)))
    new_testsuite.set('name', f"Filtered Tests - {datetime.now().isoformat()}")
    
    # Write to output file
    tree = ET.ElementTree(new_root)
    ET.indent(tree, space="  ", level=0)
    tree.write(output_file, encoding='utf-8', xml_declaration=True)
    
    # Count passed tests in each file
    base_passed = sum(1 for key in base_cases if get_test_status(base_elements[key]) == 'PASSED')
    compare_passed = sum(1 for key in compare_cases if get_test_status(compare_elements[key]) == 'PASSED')
    
    print(f"Base tests: {len(base_cases)} ({base_passed} passed), Compare tests: {len(compare_cases)} ({compare_passed} passed)")
    print(f"Created filtered JUnit XML with {len(new_or_changed_tests)} new or changed tests")
    print(f"Output saved to: {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python junit_regression.py <base_report.xml> <compare_report.xml> [output_filtered.xml]")
        sys.exit(1)
    
    base_file = sys.argv[1]
    compare_file = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else "filtered_junit.xml"
    
    create_filtered_junit_xml(base_file, compare_file, output_file)
