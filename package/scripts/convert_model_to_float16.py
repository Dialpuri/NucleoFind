import argparse
import onnx
from onnxconverter_common import float16

def convert_to_float16(input_path, output_path):
    """Converts an ONNX model to float16 format."""
    try:
        print(f"Loading model from {input_path}...")
        model = onnx.load(input_path)

        print("Converting model to float16...")
        model_fp16 = float16.convert_float_to_float16(model)

        print(f"Saving converted model to {output_path}...")
        onnx.save(model_fp16, output_path)
        print("Conversion completed successfully.")
    except Exception as e:
        print(f"An error occurred during conversion: {e}")

def main():
    parser = argparse.ArgumentParser(description="Convert an ONNX model to float16 format.")
    parser.add_argument("input", type=str, help="Path to the input ONNX model file.")
    parser.add_argument("output", type=str, help="Path to save the converted ONNX model file.")

    args = parser.parse_args()
    convert_to_float16(args.input, args.output)

if __name__ == "__main__":
    main()

