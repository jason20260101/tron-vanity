#!/bin/bash

# TRON Vanity Address Generator - Build Script
# 用于打包当前平台的发布版本

set -e

VERSION=${1:-"v1.0.0"}
PLATFORM=$(uname -s | tr '[:upper:]' '[:lower:]')
ARCH=$(uname -m)

# 确定平台名称
case "$PLATFORM" in
    darwin)
        if [ "$ARCH" = "arm64" ]; then
            PLATFORM_NAME="macos-arm64"
        else
            PLATFORM_NAME="macos-x64"
        fi
        ;;
    linux)
        PLATFORM_NAME="linux-x64"
        ;;
    mingw*|msys*|cygwin*)
        PLATFORM_NAME="windows-x64"
        ;;
    *)
        echo "未知平台: $PLATFORM"
        exit 1
        ;;
esac

echo "================================================"
echo "TRON Vanity Generator - 构建脚本"
echo "================================================"
echo "版本: $VERSION"
echo "平台: $PLATFORM_NAME"
echo ""

# 清理并编译
echo "正在编译..."
make clean
make

# 创建发布目录
RELEASE_DIR="release/tron-vanity-${VERSION}-${PLATFORM_NAME}"
rm -rf "$RELEASE_DIR"
mkdir -p "$RELEASE_DIR"

# 复制文件
echo "正在打包..."
cp tron-vanity "$RELEASE_DIR/"
cp README.md "$RELEASE_DIR/"
cp profanity.cl "$RELEASE_DIR/"
cp keccak.cl "$RELEASE_DIR/"

# 创建压缩包
cd release
if [ "$PLATFORM" = "darwin" ] || [ "$PLATFORM" = "linux" ]; then
    tar -czvf "tron-vanity-${VERSION}-${PLATFORM_NAME}.tar.gz" "tron-vanity-${VERSION}-${PLATFORM_NAME}"
else
    zip -r "tron-vanity-${VERSION}-${PLATFORM_NAME}.zip" "tron-vanity-${VERSION}-${PLATFORM_NAME}"
fi

echo ""
echo "================================================"
echo "构建完成!"
echo "输出文件: release/tron-vanity-${VERSION}-${PLATFORM_NAME}.tar.gz"
echo "================================================"

