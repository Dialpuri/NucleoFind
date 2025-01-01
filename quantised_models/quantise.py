from pathlib import Path
import onnx
from mpmath import convert
# from onnxruntime.quantization import quantize_dynamic, QuantType, quant_pre_process
from onnxconverter_common import float16


def quantize_():
    model_dir = Path("/Users/dialpuri/Development/nucleofind/models")
    output_dir = Path("/Users/dialpuri/Development/nucleofind/quantised_models")
    for model in model_dir.glob("*"):
        preprocessed_output_model = output_dir / f"{model.stem}.onnx"
        output_model = output_dir / model.name

        quant_pre_process(model, preprocessed_output_model)
        quantized_model = quantize_dynamic(preprocessed_output_model, output_model)


def convert_to_float16():
    for model_path in model_dir.glob("*"):
        output_model = output_dir / model_path.name

        model = onnx.load(model_path)
        model_fp16 = float16.convert_float_to_float16(model, min_positive_val=1e-10, max_finite_val=2e4)
        onnx.save(model_fp16, output_model)



if __name__ == "__main__":
    model_dir = Path("/Users/dialpuri/Development/nucleofind/models")
    output_dir = Path("/Users/dialpuri/Development/nucleofind/quantised_models")

    convert_to_float16()