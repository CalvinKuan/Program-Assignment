# ==========================================
#  FPGA Placement SA Algorithm Makefile
# ==========================================

# 編譯器設定
CXX      = g++
# 編譯參數
CXXFLAGS = -std=c++17 -O3 -Wall -Wextra -pedantic -DNDEBUG

# 檔案設定
TARGET   = placer
SRC      = SA_FPGA.cpp

# 路徑設定 (相對於 src 目錄)
IN_DIR   = ../testcases
OUT_DIR  = ../output

# ==========================================
#  主要規則
# ==========================================

# 預設目標：只打 make 時執行編譯
all: $(TARGET)

# 編譯執行檔
$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $^ -o $@

# 清除編譯產物
clean:
	rm -f $(TARGET)

# ==========================================
#  一次執行所有測試 (Run All)
# ==========================================
# 這會依序執行下方定義的所有 run_xxx 指令
# 若其中一個失敗，Make 會停止執行後續指令
run_all: run_alu4 run_apex4 run_elliptic run_clma run_pdc run_pdc20k
	@echo "----------------------------------------"
	@echo "All test cases finished successfully!"
	@echo "----------------------------------------"

# ==========================================
#  個別測試指令
# ==========================================

# 1. alu4
run_alu4: $(TARGET)
	@echo "Running alu4..."
	./$(TARGET) $(IN_DIR)/alu4.txt $(OUT_DIR)/alu4_output.txt

# 2. apex4
run_apex4: $(TARGET)
	@echo "Running apex4..."
	./$(TARGET) $(IN_DIR)/apex4.txt $(OUT_DIR)/apex4_output.txt

# 3. elliptic
run_elliptic: $(TARGET)
	@echo "Running elliptic..."
	./$(TARGET) $(IN_DIR)/elliptic.txt $(OUT_DIR)/elliptic_output.txt

# 4. clma (使用 clma_2.txt)
run_clma: $(TARGET)
	@echo "Running clma..."
	./$(TARGET) $(IN_DIR)/clma_2.txt $(OUT_DIR)/clma_output.txt

# 5. pdc
run_pdc: $(TARGET)
	@echo "Running pdc..."
	./$(TARGET) $(IN_DIR)/pdc.txt $(OUT_DIR)/pdc_output.txt

# 6. pdc_20000
run_pdc20k: $(TARGET)
	@echo "Running pdc_20000..."
	./$(TARGET) $(IN_DIR)/pdc_20000.txt $(OUT_DIR)/pdc_20000_output.txt


# ==========================================
#  驗證規則 (Checker)
# ==========================================

# Checker 執行檔路徑
CHECKER_EXE = ./checker

# 1. Check alu4
check_alu4:
	$(CHECKER_EXE) $(IN_DIR)/alu4.txt $(OUT_DIR)/alu4_output.txt

# 2. Check apex4
check_apex4:
	$(CHECKER_EXE) $(IN_DIR)/apex4.txt $(OUT_DIR)/apex4_output.txt

# 3. Check elliptic
check_elliptic:
	$(CHECKER_EXE) $(IN_DIR)/elliptic.txt $(OUT_DIR)/elliptic_output.txt

# 4. Check clma
check_clma:
	$(CHECKER_EXE) $(IN_DIR)/clma_2.txt $(OUT_DIR)/clma_output.txt

# 5. Check pdc
check_pdc:
	$(CHECKER_EXE) $(IN_DIR)/pdc.txt $(OUT_DIR)/pdc_output.txt

# 6. Check pdc_20000
check_pdc20k:
	$(CHECKER_EXE) $(IN_DIR)/pdc_20000.txt $(OUT_DIR)/pdc_20000_output.txt

# 一鍵執行所有 Check
check_all: check_alu4 check_apex4 check_elliptic check_clma check_pdc check_pdc20k