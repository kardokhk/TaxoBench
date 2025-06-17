#!/usr/bin/env python3

import sys
import os
sys.path.append('scripts/src')
from evaluate_results import evaluate_classification

def test_all_classifiers():
    """
    Test classification logic for all 6 classifiers in the framework:
    1. taxonova
    2. vsearch  
    3. blast
    4. naive-bayes
    5. idtaxa
    6. spingo
    """
    
    print("=" * 80)
    print("COMPREHENSIVE CLASSIFICATION TEST FOR ALL 6 CLASSIFIERS")
    print("=" * 80)
    print()
    
    # Define test cases for each classifier based on their typical output formats
    classifiers = {
        'taxonova': {
            'description': 'TaxoNova - Modern ML-based classifier',
            'format': 'Standard QIIME2 format with prefixes',
            'test_cases': [
                {
                    'name': 'Perfect match',
                    'observed': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella; s__Prevotella intermedia',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella; s__Prevotella intermedia',
                    'expected_result': 'match'
                },
                {
                    'name': 'Underclassification - missing species',
                    'observed': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella; s__Prevotella intermedia',
                    'expected_result': 'underclassification'
                },
                {
                    'name': 'Overclassification - extra specificity',
                    'observed': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella; s__Prevotella intermedia',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella',
                    'expected_result': 'overclassification'
                },
                {
                    'name': 'Misclassification - wrong genus',
                    'observed': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides fragilis',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella; s__Prevotella intermedia',
                    'expected_result': 'misclassification'
                }
            ]
        },
        
        'vsearch': {
            'description': 'VSEARCH - Sequence similarity search',
            'format': 'Standard QIIME2 format with prefixes',
            'test_cases': [
                {
                    'name': 'Perfect match',
                    'observed': 'd__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus agalactiae',
                    'expected': 'd__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus agalactiae',
                    'expected_result': 'match'
                },
                {
                    'name': 'Underclassification - Unassigned',
                    'observed': 'Unassigned',
                    'expected': 'd__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus; s__Streptococcus agalactiae',
                    'expected_result': 'underclassification'
                },
                {
                    'name': 'Genus-level match',
                    'observed': 'd__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus',
                    'expected': 'd__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; f__Streptococcaceae; g__Streptococcus',
                    'expected_result': 'match'
                }
            ]
        },
        
        'blast': {
            'description': 'BLAST - Basic Local Alignment Search Tool',
            'format': 'Standard QIIME2 format with prefixes',
            'test_cases': [
                {
                    'name': 'Perfect match',
                    'observed': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia coli',
                    'expected': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia coli',
                    'expected_result': 'match'
                },
                {
                    'name': 'No blast hit',
                    'observed': 'No blast hit',
                    'expected': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia coli',
                    'expected_result': 'underclassification'
                },
                {
                    'name': 'Wrong family',
                    'observed': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Vibrionaceae; g__Vibrio; s__Vibrio cholerae',
                    'expected': 'd__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia coli',
                    'expected_result': 'misclassification'
                }
            ]
        },
        
        'naive-bayes': {
            'description': 'Naive Bayes - Probabilistic classifier',
            'format': 'Standard QIIME2 format with prefixes',
            'test_cases': [
                {
                    'name': 'Perfect match',
                    'observed': 'd__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium; s__Bifidobacterium longum',
                    'expected': 'd__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium; s__Bifidobacterium longum',
                    'expected_result': 'match'
                },
                {
                    'name': 'Family-level classification',
                    'observed': 'd__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae',
                    'expected': 'd__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium; s__Bifidobacterium longum',
                    'expected_result': 'underclassification'
                },
                {
                    'name': 'Wrong order',
                    'observed': 'd__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Corynebacteriales; f__Corynebacteriaceae; g__Corynebacterium',
                    'expected': 'd__Bacteria; p__Actinobacteriota; c__Actinobacteria; o__Bifidobacteriales; f__Bifidobacteriaceae; g__Bifidobacterium',
                    'expected_result': 'misclassification'
                }
            ]
        },
        
        'idtaxa': {
            'description': 'IDTAXA - Iterative Decision Tree Algorithm for Taxonomic Assignment',
            'format': 'May include unclassified_ prefixes',
            'test_cases': [
                {
                    'name': 'Perfect match',
                    'observed': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides fragilis',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides; s__Bacteroides fragilis',
                    'expected_result': 'match'
                },
                {
                    'name': 'IDTAXA unclassified format',
                    'observed': 'd__Bacteria; p__unclassified_Bacteria; c__unclassified_Bacteria; o__unclassified_Bacteria; f__unclassified_Bacteria; g__unclassified_Bacteria',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides',
                    'expected_result': 'underclassification'
                },
                {
                    'name': 'Partial unclassified',
                    'observed': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__unclassified_Bacteroidia; f__unclassified_Bacteroidia; g__unclassified_Bacteroidia',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides',
                    'expected_result': 'underclassification'
                }
            ]
        },
        
        'spingo': {
            'description': 'SPINGO - Species-level Identification of Metagenomic Amplicons',
            'format': 'Only genus and species levels, empty higher levels',
            'test_cases': [
                {
                    'name': 'Perfect genus match (genus-level eval)',
                    'observed': 'd__; p__; c__; o__; f__; g__Prevotella',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella',
                    'expected_result': 'match'
                },
                {
                    'name': 'Perfect species match (species-level eval)',
                    'observed': 'd__; p__; c__; o__; f__; g__Prevotella; s__Prevotella intermedia',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella; s__Prevotella intermedia',
                    'expected_result': 'match'
                },
                {
                    'name': 'Empty genus (underclassification)',
                    'observed': 'd__; p__; c__; o__; f__; g__',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella',
                    'expected_result': 'underclassification'
                },
                {
                    'name': 'Ambiguous genus (underclassification)',
                    'observed': 'd__; p__; c__; o__; f__; g__AMBIGUOUS',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella',
                    'expected_result': 'underclassification'
                },
                {
                    'name': 'Species when expected genus (overclassification)',
                    'observed': 'd__; p__; c__; o__; f__; g__Prevotella; s__Prevotella intermedia',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella',
                    'expected_result': 'overclassification'
                },
                {
                    'name': 'Wrong genus (misclassification)',
                    'observed': 'd__; p__; c__; o__; f__; g__Bacteroides',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella',
                    'expected_result': 'misclassification'
                },
                {
                    'name': 'Correct genus, wrong species (misclassification)',
                    'observed': 'd__; p__; c__; o__; f__; g__Prevotella; s__Prevotella melaninogenica',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Prevotella; s__Prevotella intermedia',
                    'expected_result': 'misclassification'
                },
                {
                    'name': 'Complex species name',
                    'observed': 'd__; p__; c__; o__; f__; g__Bacteroides_H_857956; s__Bacteroides_H_857956 uniformis',
                    'expected': 'd__Bacteria; p__Bacteroidota; c__Bacteroidia; o__Bacteroidales; f__Bacteroidaceae; g__Bacteroides_H_857956; s__Bacteroides_H_857956 uniformis',
                    'expected_result': 'match'
                }
            ]
        }
    }
    
    # Run tests for each classifier
    total_tests = 0
    total_passed = 0
    total_failed = 0
    failed_tests = []
    
    for classifier_name, classifier_info in classifiers.items():
        print(f"🔍 TESTING {classifier_name.upper()}")
        print(f"Description: {classifier_info['description']}")
        print(f"Format: {classifier_info['format']}")
        print("-" * 60)
        
        classifier_passed = 0
        classifier_failed = 0
        
        for i, test_case in enumerate(classifier_info['test_cases'], 1):
            total_tests += 1
            
            # Run the classification
            actual_result = evaluate_classification(
                test_case['observed'], 
                test_case['expected']
            )
            
            # Check if result matches expectation
            if actual_result == test_case['expected_result']:
                status = "✅ PASS"
                classifier_passed += 1
                total_passed += 1
            else:
                status = "❌ FAIL"
                classifier_failed += 1
                total_failed += 1
                failed_tests.append({
                    'classifier': classifier_name,
                    'test': test_case['name'],
                    'expected': test_case['expected_result'],
                    'actual': actual_result,
                    'observed': test_case['observed'],
                    'expected_tax': test_case['expected']
                })
            
            print(f"{i:2d}. {test_case['name']:<40} {status}")
            if actual_result != test_case['expected_result']:
                print(f"    Expected: {test_case['expected_result']:<20} Actual: {actual_result}")
        
        print(f"\n{classifier_name.upper()} Results: {classifier_passed} passed, {classifier_failed} failed")
        print("=" * 80)
        print()
    
    # Print summary
    print("FINAL SUMMARY")
    print("=" * 80)
    print(f"Total tests: {total_tests}")
    print(f"Passed: {total_passed}")
    print(f"Failed: {total_failed}")
    print(f"Success rate: {(total_passed/total_tests)*100:.1f}%")
    
    if failed_tests:
        print("\n❌ FAILED TESTS:")
        print("-" * 40)
        for i, failed in enumerate(failed_tests, 1):
            print(f"{i}. {failed['classifier']} - {failed['test']}")
            print(f"   Expected: {failed['expected']}")
            print(f"   Actual: {failed['actual']}")
            print(f"   Observed: {failed['observed']}")
            print(f"   Expected taxonomy: {failed['expected_tax']}")
            print()
    else:
        print("\n🎉 ALL TESTS PASSED!")
        print("Classification logic is consistent across all 6 classifiers!")
    
    print("=" * 80)

if __name__ == "__main__":
    test_all_classifiers()
