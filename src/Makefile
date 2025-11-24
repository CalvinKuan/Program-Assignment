# 編譯器
CXX      = g++
# 編譯參數（如果要方便 debug 可以把 -DNDEBUG 拿掉，改成 -g）
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -DNDEBUG

# 目標執行檔名稱
TARGET   = placer
# 原始碼
SRC      = SA_FPGA.cpp

# 預設目標：編譯出執行檔
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

# 方便測試用：make run IN=input.txt OUT=out.txt
run1: $(TARGET)
	./$(TARGET) ../testcases/alu4.txt ../output/alu4output.txt
run5:
	./$(TARGET) ../testcases/pdc.txt ../output/pdcoutput.txt

# 清除編譯產物
clean:
	rm -f $(TARGET)
